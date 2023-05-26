/************************************************************************
 * RIA, version 1.1
 * Copyright 2015-2023
 * Richard Howey
 * Institute of Genetic Medicine, Newcastle University
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


/*! \file Data.h
    \brief This file contains classes for manipulating SNP data.
    
*/

#ifndef __DATA
#define __DATA

#include <set>
#include <map>
#include <list>
#include <stack>
#include <fstream>
#include <string>

using namespace std;

class SNPWindow;
class Analysis;

//! Stores SNP data for all subjects for a given SNP.
struct SNPData
{
	
	list<unsigned char> noMinorAllelesAllSubjects; //ordered list of the number of minor alleles for the subjects	
	string name;
	double cM;
	unsigned int bpPos;
	unsigned int chromosome;

	SNPData(SNPWindow * snpWindow, const string & nm, const double & cm, const unsigned int & bp, const unsigned int & chr);	

	~SNPData()
	{


	};

	void readInNewData(SNPWindow * snpWindow, const string & nm, const double & cm, const unsigned int & bp, const unsigned int & chr);
};

//! Stores the IBD estimates. 
struct IBDData
{
	double f1j; //est. prob of sharing 1 allele IBD, f0j = 1 - f1j -f2j
	double f2j; //est. prob of sharing 2 alleles IBD
	bool error;

	IBDData(const double & f1, const double & f2, const bool & er) : f1j(f1), f2j(f2), error(er) {};

	~IBDData() {};
};

//! Stores IBD estimates for all affected relative pairs (ARPs), could be priors or posteriors.
struct IBDEstimates
{
	map<string, IBDData *> snpIBDEstimates; //string = famID+indivID1+indivID2, IBD data estimated probs of alleles IBD
	
	IBDEstimates() : snpIBDEstimates() {};
	
	~IBDEstimates()
	{
		for(map<string, IBDData *>::iterator i = snpIBDEstimates.begin(); i != snpIBDEstimates.end(); ++i)
		{
			delete i->second;
		};
	};
};

//! Stores SNP data for the SNP window in question.
class SNPWindow
{
private:

	list<double> geneticDistances; //genetic Distances of all SNPs in order
	list<unsigned int> basePairs; //base-pair positions of all SNPs in order
	list<unsigned int> chromosomes; //chromosomes of all SNPs in order
	list<string> namesOfSNPs;	
	map<unsigned int, SNPData *> window; //SNP no, SNP Data
	stack<SNPData *> spareSNPs;
	double windowCMSize; //window size for full window cM
	bool missingCMZero;
	bool decreasingCM;
	unsigned int windowStepSize;
	unsigned int totalNoSubjects;
	list<unsigned char> caseControls; //list of bool values to whether a subject is a case or not	
	bool newChr;
	string cpStr, rmStr, endCommand, postFamFile;
	unsigned int startSNP;
	unsigned int endSNP;
	unsigned int startSNPJobNo;
	unsigned int endSNPJobNo;
	double noCalcs;

	SNPData * analysisSNP;

	unsigned int analysisSNPNo;
	unsigned int windowMinSNPSize;
	bool validWindow;
	double firstCMChr;

	list<unsigned int> windowSizes;
	list<double>::const_iterator analysisCMSNP; //will be in middle of the SNP window
	list<unsigned int>::const_iterator analysisChrSNP; //will be in middle of the SNP window
	list<double>::const_iterator endCMSNP; //the endSNP iterator is the SNP after the SNP window (i.e. is not in window)
	list<unsigned int>::const_iterator endChr; //the endChr iterator tracks the chr of the SNP window end, needed for chr switching
	list<string>::const_iterator nameSNPiter; //used for adding new SNP data
	list<double>::const_iterator cmSNPiter; //used for adding new SNP data
	list<unsigned int>::const_iterator bpSNPiter; //used for adding new SNP data
	list<unsigned int>::const_iterator chrSNPiter; //used for adding new SNP data
	unsigned int addingSNPNo; //used for adding new SNP data
	bool newSNPObject;	

	//for reading in data from binary files
	ifstream readSNPData;
	unsigned int bitCount;
	int one;
	int aBit;	
	char buffer[1];

	Analysis * analysis;

public:

	SNPWindow(string & filename, double & wcMs, unsigned int & wmssnps, unsigned int & wstps, bool & misscm, bool decm, unsigned int & stSNPno, unsigned int & eSNPno, string & stSNPnoname, string & eSNPname, string & cp, string & rm, string & end, string & pff, const bool & priorsGiven, bool & priorOnly, unsigned int & jobNo, unsigned int & jobTotal, Analysis * analysis);

	~SNPWindow()
	{
			
		readSNPData.close();		

		while(!spareSNPs.empty())
		{
			delete spareSNPs.top();
			spareSNPs.pop();
		};
	};
	
	//methods regarding the SNP descriptions
	void setUpSNPDesciptionData(string & filename, string & startSNPName, string & endSNPName);
	void updateStartAndEndSNP(unsigned int & startSNPAnal, unsigned int & endSNPAnal);
	void setWindowMinSNPSize(unsigned int & wmss) {windowMinSNPSize = wmss;};
	double getWindowCMSize() {return windowCMSize;};
	void setPostFamFilename(string & cp, int & pid);
	bool updateEndSNP(); 
	void updateStartSNP();
	unsigned int getNoSNPs() {return geneticDistances.size();};
	void displayWindowStats();
	void recordWindowSize() {windowSizes.push_back(window.size());};

	void setStartAndEndSNPBasedOnJobNo(unsigned int & jobNo, unsigned int & jobTotal);
	bool isAnalysisSNPInJobWindow();
	void setUpFirstWindow(string & filename);
	void setUpStartAndEndSNPs();
	bool nextWindow();	
	void addNextSNPDataToWindow();
	bool isWindowInOneChr();
	bool hasWindowMinNoSNPs() {return window.size() >= windowMinSNPSize;};
	unsigned int getLastWindowChr();
	bool getValidWindow() {return validWindow;};
	void setValidWindow(const bool & vw) {validWindow = vw;};

	SNPData * getAnalysisSNPData() const {return analysisSNP;};

	unsigned char getNextNoOfMinorAlleles();
	void setupCaseControls(string & fname, const bool & makePostFamFile);
	unsigned int getTotalNoSubjects() {return totalNoSubjects;};
	unsigned int getStartCMWindowSNP(unsigned int & startAnalysisSNPno);
	void advanceToFirstWindow(unsigned int & analysisSNPNo);
	void advanceToFirstWindowInChr();
	void advanceSNPData();
	SNPData * getNewSNPData();
	
	double getAnalysisSNPCM() const {return analysisSNP->cM;};
	unsigned int getAnalysisSNPBP() const {return analysisSNP->bpPos;};
	unsigned int getAnalysisChromosome() const {return analysisSNP->chromosome;};
	string getAnalysisSNPName() const {return analysisSNP->name;};
	unsigned int getAnalysisSNPNo() const {return analysisSNPNo;};
	void deleteAllTemporaryFiles();
	void displayWindowInfo();
	void outputProgress(const unsigned int & calcCount);
	void estimateTotalPostCalcs(const unsigned int & jobTotal);

	void createFilesForPosteriorCalc(string & outputFilename);
	void createFamFileForPosteriorCalc(string & outputFilename);
	void createBimFileForPosteriorCalc(string & outputFilename);
};

#endif
