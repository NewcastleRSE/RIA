/************************************************************************
 * RIA, version 1.1
 * Copyright 2015-present
 * Richard Howey
 * Research Software Engineering, Newcastle University
 *
 * richard.howey@ncl.ac.uk
 * http://www.staff.ncl.ac.uk/richard.howey/
 *
 * This file is part of RIA, the pedigree file processing program for EMIM.
 *
 * RIA is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * RIA is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with RIA.  If not, see <http://www.gnu.org/licenses/>.
 ************************************************************************/


/*! \file Analysis.h
    \brief This file organises the analysis.
    
*/

#ifndef __ANALYSIS
#define __ANALYSIS

#include <string>
#include <iostream>
#include <sstream>

using namespace std; // initiates the "std" or "standard" namespace

#include "main.h"
#include "Data.h"
#include "Model.h"
#include "ModelFitting.h"

#ifdef USING_WINDOWS
#include "process.h"
#endif

#ifndef USING_WINDOWS
#include <sys/types.h>
#include <unistd.h>
#endif

//! Returns a string of the run time.
string getTime(const double & t);

//! Class for storing three parameters.
struct threeParameters
{
	double beta0, beta1, beta2;

	threeParameters() : beta0(0), beta1(0), beta2(0) {};
	threeParameters(const double & v1, const double & v2, const double & v3) : beta0(v1), beta1(v2), beta2(v3) {};
};

//! Organises the correect analysis to perform.
class Analysis
{
private:

	//parameters for the options of the analysis to perform
	string filename;
	string outputFilename;
	string plink;
	string king;
	string truffle;
	string trufflePairsFilename;
	string truffleOptions;
	string gzip;
	string ibdStitch;
	string rmStr, cpStr, endCommand;
	int pid;
	string postFamFile;
	
	string posteriorInputFilePrefix;
	string posteriorOutputFilePrefix;
	unsigned int posteriorStartWindow;
	unsigned int posteriorEndWindow;

	unsigned int windowStepSize;
	bool missingCMZero;
	bool decreasingCM;
	unsigned int startSNP;
	unsigned int endSNP;
	unsigned int jobNo;
	unsigned int jobTotal;
	
	string priorFilename;
	string priorOutputFilename;
	string plinkPriorOptions;
	IBDEstimates * priors; 
	IBDEstimates * posteriors;

	//column numbers
	unsigned int famIDCol, indivID1Col, indivID2Col, noSNPsCol, z0Col, phiCol, IBD0Col, kinshipCol, errorCodeCol;
	unsigned int noColumns;

	bool priorOnly;
	bool noDomVar;
	bool doIbdStitch;
	bool doTruffle;
	bool keepTempFiles;
	double prevVarA;
	double prevVarD;

	string outputPriorFilename1;
	string outputPriorFilename2;
	string outputPriorFilename3IBDStitch;
	string outputPriorFilename3Truffle;
	string outputPostFilename;

	SNPWindow * snpWindow;
	ModelIBD * modelIBD;
	
	ofstream resultsFile;
	
	//used for reading in data
	char oneBuffer[1];
	int aBit;

	static const double oneOverSqRoot2;

public:

	Analysis(string & fn, string & ofn, string & gfn, string & gofn, string & plk, string & plkops, string & kng, string & trff, string & trffOps, string & gzp, string & ibdSt, double & wcMs, bool & misscm, bool decm, unsigned int & mwssnps, unsigned int & wstps, unsigned int & ssnp, unsigned int & esnp, string & ssnpname, string & esnpname, unsigned int & jn, unsigned int & jt, bool & ndv, bool & go, bool & ktf,
		string & pifp, string & pofp, unsigned int & psw, unsigned int & pew)
		:  filename(fn), outputFilename(ofn), priorFilename(gfn), priorOutputFilename(gofn), plink(plk), plinkPriorOptions(plkops), king(kng), truffle(trff), truffleOptions(trffOps), gzip(gzp), ibdStitch(ibdSt), windowStepSize(wstps), missingCMZero(misscm), decreasingCM(decm), startSNP(ssnp), endSNP(esnp), jobNo(jn), jobTotal(jt), noDomVar(ndv), priorOnly(go), prevVarA(0), prevVarD(0), keepTempFiles(ktf),
		posteriorInputFilePrefix(pifp), posteriorOutputFilePrefix(pofp), posteriorStartWindow(psw), posteriorEndWindow(pew)
	{			
		#ifndef USING_WINDOWS
			pid = getpid();
			endCommand = " >/dev/null 2>&1";
			rmStr = "rm ";
			cpStr = "cp ";
		#endif

		#ifdef USING_WINDOWS
			pid = _getpid();
			endCommand = " >nul 2>&1";
			rmStr = "del ";
			cpStr = "copy ";
		#endif

		doTruffle = truffle != "";
		doIbdStitch = ibdStitch != "";
		setupTemporaryFilenames();

		//setup priors and posteriors
		priors = new IBDEstimates();
		posteriors = new IBDEstimates();

		//setup the model
		modelIBD = new ModelIBD();

		//add pointers in model to priors and posteriors
		modelIBD->setPriorsAndPosteriors(priors, posteriors);

		//setup SNP window
		snpWindow = new SNPWindow(filename, wcMs, mwssnps, windowStepSize, missingCMZero, decreasingCM, startSNP, endSNP, ssnpname, esnpname, cpStr, rmStr, endCommand, postFamFile, (priorFilename != ""), priorOnly, jobNo, jobTotal, this); 
		resultsFile.unsetf(ios::floatfield);            
		resultsFile.precision(10);

		snpWindow->updateStartAndEndSNP(startSNP, endSNP); //ensure SNP numbers match if updated with the SNP name

		if(doTruffle) outputTrufflePairsFile();

		if(priorFilename != "") setupPriorsUsingInputFile(); //read in priors from file
		else calculatePriors(); //calculate priors using prior estimate of IBDs

	};

	//! Delete analysis things
	~Analysis()
	{		
		if(!priorOnly) delete snpWindow;
		delete modelIBD;
		delete priors;
		delete posteriors;
	};

	void runAnalysis();	
	void pruneDataForPriorCalc();
	void calculatePriors();
	void setupPriors(string & filenameIBD);
	void setupPriorsUsingInputFile();
	void setupPriorsTruffle(string & filenameIBD);
	void outputTrufflePairsFile();
	void calculatePosteriors(unsigned int & windowNumber);
	void setupPosteriors(string & filenameIBD);
	void setupPosteriorsTruffle(string & filenameIBD);
	void fitModels(const unsigned int & snpID);
	void setupTemporaryFilenames();
	void deletePriorTemporaryFiles();
	void deletePostTemporaryFiles();
	void deletePostFamTemporaryFile();
	void deletePairsTemporaryFile();
	void deleteAllTemporaryFiles();
	void checkKingKinFile(string & filenameIBD);

	void calculatePriorsTruffle();
	void calculatePosteriorsTruffle();

	void calculatePriorsIBDStitch();
	void calculatePosteriorsIBDStitch();

	unsigned int getNextNoOfMinorAlleles(ifstream & readSNPData, unsigned int & bitCount);
	double getPvalueChiSq1DF(double & chisq);
};


#endif
