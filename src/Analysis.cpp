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


/*! \file Analysis.cpp
    \brief This file contains the methods the various analyse.
    
*/
#include <iostream>
#include <sstream>
#include <math.h>
#include <cstdlib>

using namespace std; // initiates the "std" or "standard" namespace

#include "Analysis.h"
#include "Model.h"
#include "Data.h"
#include "ModelFitting.h"
#include "main.h"
#include "cdflib.h"


#ifdef USING_WINDOWS
	#include <process.h>
#endif

const double Analysis::oneOverSqRoot2 = 0.7071067811865475244008443621048490392848359376884740365883;

//! Converts an integer to a string.
string toString(int & i)
{
	ostringstream aStringStream;
	aStringStream << i;

	return aStringStream.str();
};

//! Converts an integer to a string.
string toStringUI(unsigned int & i)
{
	ostringstream aStringStream;
	aStringStream << i;

	return aStringStream.str();
};

//! Returns a string of the run time.
string getTime(const double & t)
{
	double time = t;
	int days = 0;
	int hours = 0;
	int minutes = 0;
	int seconds = 0;

	string ans = "";
	days = (int)(time / 86400); time -= days*864000.0;

	hours = (int) (time / 3600); time -= hours*3600;
	minutes = (int) (time / 60); time -= minutes*60;
	seconds = (int) time;

	if(days == 1) ans += "1 day";
	else if(days > 0) ans += toString(days) + " days";

	if(hours > 0)
	{
		if(days != 0)
		{
			if(minutes == 0 && seconds == 0) ans += " and ";
			else ans += ", ";
		};

		if(hours == 1) ans += "1 hour";
		else ans += toString(hours) + " hours";
	};

	if(minutes > 0)
	{
		if(ans != "")
		{
			if(seconds == 0) ans += " and ";
			else ans += ", ";
		};

		if(minutes == 1) ans += "1 minute";
		else ans += toString(minutes) + " minutes";
	};

	if(seconds > 0)
	{
		if(ans != "")
		{
			ans += " and ";			
		};

		if(seconds == 1) ans += "1 second";
		else ans += toString(seconds) + " seconds";
	};

	if(ans == "") ans = "less than one second";

	return ans;
};

//! Gets the log filename from the results file name.
string getSNPNotAnalFileName(string & outFileName)
{
	unsigned int filenameLength = outFileName.length();
	string notAnalFileName;

	//find extension
	unsigned int a = filenameLength - 1;
	while(a > 0)
	{
		if(outFileName.substr(a, 1) == ".") break;
		a--;
	};

	if(a > 0) notAnalFileName = outFileName.substr(0, a) + "-SNPs-not-anal.txt";
	else notAnalFileName = outFileName + "-SNPs-not-anal.txt";

	return notAnalFileName;
};

//! Sets up files names used for temporary files.
void Analysis::setupTemporaryFilenames()
{
	//convert PID to string, use PID to create unique temporary file names in case many progs are ran
	ostringstream aStringStream;
	aStringStream << pid;
	string pidStr = aStringStream.str();

	outputPriorFilename1 = "tempRIA-priors1-"+pidStr;
	outputPriorFilename2 = "tempRIA-priors2-"+pidStr;
	outputPriorFilename3IBDStitch = "tempRIA-priors3-" + pidStr;
	outputPriorFilename3Truffle = "tempRIA-priors3-" + pidStr;

	outputPostFilename = "tempRIA-posterior-" + pidStr;

	trufflePairsFilename = "tempRIA-pairs-" + pidStr + ".dat";

	//setup fam file name needed for posterior calcs
	postFamFile = "tempRIA-post-"+pidStr + ".fam"; //the same as set in Analysis;

};

//! Deletes all temporary files.
void Analysis::deleteAllTemporaryFiles()
{
	deletePriorTemporaryFiles();
	deletePostTemporaryFiles();
	deletePostFamTemporaryFile();
	if(doTruffle) deletePairsTemporaryFile();
};

//! Deletes prior files used for temporary files.
void Analysis::deletePriorTemporaryFiles()
{
	string rmStrPri1 = rmStr + outputPriorFilename1 + "*" + endCommand;
	string rmStrPri2 = rmStr + outputPriorFilename2 + "*" + endCommand;
	
	//delete temp files
	if(!keepTempFiles) system(rmStrPri1.c_str());
	if(!keepTempFiles) system(rmStrPri2.c_str());
	if (doTruffle && !keepTempFiles)
	{
		string rmStrPri3 = rmStr + outputPriorFilename3Truffle + "*" + endCommand;
		system(rmStrPri3.c_str());		
	}
	else if(doIbdStitch && !keepTempFiles)
	{
		string rmStrPri3 = rmStr + outputPriorFilename3IBDStitch + "*" + endCommand;
		system(rmStrPri3.c_str());
		string rmStrPri4 = rmStr + "tempRIA-paras" + toString(pid) + "*" + endCommand;
		system(rmStrPri4.c_str());
	};
};

//! Deletes posterior files used for temporary files.
void Analysis::deletePostTemporaryFiles()
{
	if(posteriorInputFilePrefix != "") return; //nothing to delete if inputting posteriors

	string rmStrPost = rmStr + outputPostFilename + "*" + endCommand;
	
	//do not want to delete *.kin posterior files if we are outputting them
	if(posteriorOutputFilePrefix != "") rmStrPost = rmStr + outputPostFilename + ".bed " + outputPostFilename + ".bim " + outputPostFilename + ".fam " + outputPostFilename + ".kin0" + endCommand;

	//delete temp files
	if(!keepTempFiles) system(rmStrPost.c_str());
};

//! Deletes posterior files used for temporary files.
void Analysis::deletePostFamTemporaryFile()
{
	//delete temp .fam file
	string rmFamFile = rmStr + postFamFile + endCommand;
	if(!keepTempFiles) system(rmFamFile.c_str());
};

//! Deletes pairs file used by truffle.
void Analysis::deletePairsTemporaryFile()
{
	//delete pairs file
	string rmPairsFile = rmStr + trufflePairsFilename + endCommand;
	if(!keepTempFiles) system(rmPairsFile.c_str());
};

//! Runs the chosen analysis.
void Analysis::runAnalysis()
{

	//priors calculated within constructor
	unsigned int noSNPsNotAnalysed = 0;	
	ofstream notAnalSNPsFile;
	
	resultsFile.open(outputFilename.c_str());
	
	//output header line for results file, but not for jobs 2 or more, if jobs used
	if(jobNo <= 1) resultsFile << "SNP CHR ID CM BP VAR_A VAR_D MLS\n";

	bool atEndOfData = false;

	//loop thro' all the SNPs updating the SNP window as we go along
	unsigned int snpID = snpWindow->getAnalysisSNPNo();
	unsigned int snpIDOld = 0;
	unsigned int calcCount = 0;
	unsigned int windowNumber = 0;

	snpWindow->outputProgress(calcCount);

	snpWindow->estimateTotalPostCalcs(jobTotal);

	while(!atEndOfData)
	{
		snpWindow->outputProgress(calcCount);

		//perform analysis for the next window
		if(snpWindow->isAnalysisSNPInJobWindow())
		{
			if(snpWindow->getValidWindow() && snpWindow->hasWindowMinNoSNPs())
			{
				windowNumber++;
				if((posteriorStartWindow == 0 || windowNumber >= posteriorStartWindow) && (posteriorEndWindow == 0 || windowNumber <= posteriorEndWindow)) //do selected windows if set
				{					
					calculatePosteriors(windowNumber);
					fitModels(snpID);
					snpWindow->recordWindowSize();
				};

				//cout << "*";				
			};

			++calcCount;
		};

		//snpWindow->displayWindowInfo(); //testing

		//move to next window
		atEndOfData = snpWindow->nextWindow();
		snpID = snpWindow->getAnalysisSNPNo();

		if(snpID == snpIDOld)
		{
			outErr("\nProblem updating SNP window, please check centimorgan SNP data is correct!\n");
			deleteAllTemporaryFiles();
			exit(1);
		};

		snpIDOld = snpID;
	};

	screenOut("\r                             "); screenFlush("");

	resultsFile.close();
	if(noSNPsNotAnalysed != 0) notAnalSNPsFile.close();
	
	snpWindow->displayWindowStats();

	deletePostFamTemporaryFile();
	if(doTruffle) deletePairsTemporaryFile();
};

//! Prunes data using PLINK for use to calculate the prior.
void Analysis::pruneDataForPriorCalc()
{
	//firstly get list of SNPs for pruned list, filter out cases and only include common SNPs, MAF > 0.4
	string prefix = filename.substr(0, filename.length() - 4);
	string plinkCommand;

	if(plinkPriorOptions != "") plinkCommand = plink + " --noweb " + plinkPriorOptions + " --bfile " + prefix
		+ " --out " + outputPriorFilename1 + endCommand;
	else plinkCommand = plink + " --noweb --indep 50 50 2 --mind 0.01 --maf 0.25 --geno 0.05 --bfile " + prefix
		+ " --out " + outputPriorFilename1 + endCommand;

	out("Creating list of pruned SNPs using PLINK command:\n"); out(plinkCommand); out("\n\n");

	unsigned int plinkErrorCode = system(plinkCommand.c_str());

	if(plinkErrorCode != 0)
	{
		outErr("Problem executing PLINK (error code "); outErr(plinkErrorCode); outErr(") when creating data file with:\n\n ");
		outErr(plinkCommand); outErr("\n\n");
		outErr("Log file: "); outErr(outputPriorFilename1); outErr(".log:\n\n");
		string plinkLog = "more " + outputPriorFilename1 + ".log";
		system(plinkLog.c_str());

		deleteAllTemporaryFiles();
		exit(plinkErrorCode);
	};

	//now get a .bed file of the pruned SNPs found in the last step
	string plinkCommand2 = plink + " --noweb --bfile " + prefix + " --filter-cases --extract " + outputPriorFilename1 + ".prune.in --make-bed --out " + outputPriorFilename2 + endCommand;

	out("Creating data file to calculate priors using PLINK command:\n"); out(plinkCommand2); out("\n\n");

	plinkErrorCode = system(plinkCommand2.c_str());

	if(plinkErrorCode != 0)
	{
		outErr("Problem executing PLINK (error code "); outErr(plinkErrorCode); outErr(") when creating data file with:\n\n ");
		outErr(plinkCommand2); outErr("\n\n");
		outErr("Log file: "); outErr(outputPriorFilename2); outErr(".log:\n\n");
		string plinkLog = "more " + outputPriorFilename2 + ".log";
		system(plinkLog.c_str());

		deleteAllTemporaryFiles();
		exit(plinkErrorCode);
	};

};

//! Calculates the priors using global IBD calculations on a reduced SNP set.
void Analysis::calculatePriors()
{
	pruneDataForPriorCalc();

	if(doTruffle) { calculatePriorsTruffle(); return; }
	else if(doIbdStitch) { calculatePriorsIBDStitch(); return; };

	//calculate the IBDs for the pruned data set in order to calculate the priors
	string kingCommand = king + " -b " + outputPriorFilename2 + ".bed --homog --prefix " + outputPriorFilename2 + endCommand;

	out("Calculating priors using KING command:\n"); out(kingCommand); out("\n\n");

	unsigned int kingErrorCode = system(kingCommand.c_str());

	if(kingErrorCode != 0)
	{
		outErr("Problem executing KING (error code "); outErr(kingErrorCode); outErr(") when estimating IBDs with:\n\n ");
		outErr(kingCommand); outErr("\n\n");

		deleteAllTemporaryFiles();
		exit(kingErrorCode);		
	};

	//make copy of .fam file for use with posterior calculations
	if(!priorOnly)
	{
		string copyFamCommand = cpStr + outputPriorFilename2 + ".fam " + postFamFile + endCommand;
		unsigned int cpErrorCode = system(copyFamCommand.c_str());
		if(cpErrorCode != 0)
		{
			outErr("Problem copying (error code "); outErr(cpErrorCode); outErr(") .fam file for posterior calculations with:\n\n ");
			outErr(copyFamCommand); outErr("\n\n");

			deleteAllTemporaryFiles();
			exit(cpErrorCode);		
		};
	};

	//read in data and store priors
	setupPriors(outputPriorFilename2);

	out("Number of affected relative pairs (ARPs) in priors: "); out(priors->snpIBDEstimates.size()); out("\n\n");

	ostringstream aStringStream;
	aStringStream << pid;
	string pidStr = aStringStream.str();

	if(!priorOnly)
	{
		//just output the posterior command here once for convience, so is not outputted for every analysis SNP
		kingCommand = king + " -b " + outputPostFilename + ".bed --homog --prefix " + outputPostFilename + endCommand;

		out("Calculating posteriors (for each SNP window) using KING command:\n"); out(kingCommand); out("\n\n");
	};

	deletePriorTemporaryFiles();
};

//! Outputs Truffle pairs file.
void Analysis::outputTrufflePairsFile()
{
	
	//try and find the family file and read in data
	unsigned int length = filename.length();
	string famFilename = filename.substr(0, length-4) + ".fam";

	ifstream readFamilyFile;
	readFamilyFile.open(famFilename.c_str());
	if(!readFamilyFile.is_open())
	{
		outErr("\nCannot read family file: "); outErr(famFilename); outErr("!\n");
		deleteAllTemporaryFiles();
		exit(1);
	};


	string famID, indivID, FatherId, MotherID, sexID, famIndivID, phenoType;
	string prevFamIndivID = "";
	map<string, set<string> > families;
	map<string, set<string> >::iterator fa;
	set<string> aFam;

	//loop thro' subjects and store the cases
	do{		
		readFamilyFile >> famID >> indivID >> FatherId >> MotherID >> sexID >> phenoType;
		famIndivID = famID + "-" + indivID;

		//do not duplicate the last row
		if(famIndivID != prevFamIndivID) 
		{
			fa = families.find(famID);
			if(fa != families.end())
			{
				fa->second.insert(indivID);
			}
			else {
				aFam.clear();
				aFam.insert(indivID);
				families[famID] = aFam;
			};
		};

		prevFamIndivID = famIndivID;
	}while(!readFamilyFile.eof());

	readFamilyFile.close();

	ofstream pairsFile(trufflePairsFilename.c_str());
	unsigned int i1Num, i2Num;
	string indiv1Str, indiv2Str;

	//output every pair in each family to pairs file.
	for(map<string, set<string> >::const_iterator fa = families.begin(); fa != families.end(); ++fa)
	{
		i1Num = 1;
		for(set<string>::const_iterator i1 = fa->second.begin(); i1 != fa->second.end(); ++i1, ++i1Num)
		{
			i2Num = 1;
			for(set<string>::const_iterator i2 = fa->second.begin(); i2 != fa->second.end(); ++i2, ++i2Num)
			{
				if(i2Num > i1Num)
				{
					indiv1Str = fa->first + "_" + *i1;
					indiv2Str = fa->first + "_" + *i2;
					pairsFile << indiv1Str << " " << indiv2Str << "\n";
				};
			};
		};
	};

	pairsFile.close();
};

//! Calculates the priors using global IBD calculations on a reduced SNP set using TRUFFLE.
void Analysis::calculatePriorsTruffle()
{
	//firstly convert to a VCF gzipped file for use with truffle
	string plinkCommand3 = plink + " --noweb --bfile " + outputPriorFilename2 + " --recode vcf --out " + outputPriorFilename3Truffle + endCommand;

	out("Creating data file to calculate priors using PLINK command:\n"); out(plinkCommand3); out("\n\n");

	unsigned int plinkErrorCode = system(plinkCommand3.c_str());

	if(plinkErrorCode != 0)
	{
		outErr("Problem executing PLINK (error code "); outErr(plinkErrorCode); outErr(") when creating data file with:\n\n ");
		outErr(plinkCommand3); outErr("\n\n");
		outErr("Log file: "); outErr(outputPriorFilename3Truffle); outErr(".log:\n\n");
		string plinkLog = "more " + outputPriorFilename3Truffle + ".log";
		system(plinkLog.c_str());

		deleteAllTemporaryFiles();
		exit(plinkErrorCode);
	};

	//now gzip the vcf file
	string gzipCommand = gzip + " -f " + outputPriorFilename3Truffle + ".vcf" + endCommand;

	unsigned int gzipErrorCode = system(gzipCommand.c_str());

	if(gzipErrorCode != 0)
	{
		outErr("Problem executing gzip (error code "); outErr(gzipErrorCode); outErr(") when gzipping data file with:\n\n ");
		outErr(gzipCommand); outErr("\n\n");

		deleteAllTemporaryFiles();
		exit(gzipErrorCode);
	};

	//calculate the IBDs for the pruned data set in order to calculate the priors
	string truffleCommand = truffle + " --pairs-file " + trufflePairsFilename + " --vcf " + outputPriorFilename3Truffle + ".vcf.gz " + truffleOptions + " --out " + outputPriorFilename2 + endCommand;

	out("Calculating priors using TRUFFLE command:\n"); out(truffleCommand); out("\n\n");

	unsigned int truffleErrorCode = system(truffleCommand.c_str());

	if(truffleErrorCode != 0)
	{
		outErr("Problem executing TRUFFLE (error code "); outErr(truffleErrorCode); outErr(") when estimating IBDs with:\n\n ");
		outErr(truffleCommand); outErr("\n\n");

		deleteAllTemporaryFiles();
		exit(truffleErrorCode);
	};

	//make copy of .fam file for use with posterior calculations
	if(!priorOnly)
	{
		string copyFamCommand = cpStr + outputPriorFilename2 + ".fam " + postFamFile + endCommand;
		unsigned int cpErrorCode = system(copyFamCommand.c_str());
		if (cpErrorCode != 0)
		{
			outErr("Problem copying (error code "); outErr(cpErrorCode); outErr(") .fam file for posterior calculations with:\n\n ");
			outErr(copyFamCommand); outErr("\n\n");

			deleteAllTemporaryFiles();
			exit(cpErrorCode);
		};
	};

	//read in data and store priors
	setupPriorsTruffle(outputPriorFilename2);

	out("Number of affected relative pairs (ARPs) in priors: "); out(priors->snpIBDEstimates.size()); out("\n\n");

	ostringstream aStringStream;
	aStringStream << pid;
	string pidStr = aStringStream.str();

	if(!priorOnly)
	{
		//just output the posterior command here once for convience, so is not outputted for every analysis SNP		
		truffleCommand = truffle + " --pairs-file " + trufflePairsFilename + " --vcf " + outputPostFilename + ".vcf.gz " + truffleOptions + " --out " + outputPostFilename + " " + endCommand;

		out("Calculating posteriors (for each SNP window) using TRUFFLE command:\n"); out(truffleCommand); out("\n\n");
	};

	deletePriorTemporaryFiles();
};

//! Calculates the priors using global IBD calculations on a reduced SNP set for IBD stitch.
void Analysis::calculatePriorsIBDStitch()
{
	//make .ped file this alleles recoded as 1 and 2
	string plinkCommand3 = plink + " --noweb --bfile " + outputPriorFilename2 + " --recode12 --out " + outputPriorFilename3IBDStitch + endCommand;

	out("Creating data file to calculate priors using PLINK command:\n"); out(plinkCommand3); out("\n\n");

	unsigned int plinkErrorCode = system(plinkCommand3.c_str());

	if(plinkErrorCode != 0)
	{
		outErr("Problem executing PLINK (error code "); outErr(plinkErrorCode); outErr(") when creating data file with:\n\n ");
		outErr(plinkCommand3); outErr("\n\n");
		outErr("Log file: "); outErr(outputPriorFilename3IBDStitch); outErr(".log:\n\n");
		string plinkLog = "more " + outputPriorFilename3IBDStitch + ".log";
		system(plinkLog.c_str());

		deleteAllTemporaryFiles();
		exit(plinkErrorCode);
	};

	//make .frq file with plink
	string plinkCommand4 = plink + " --noweb --file " + outputPriorFilename3IBDStitch + " --freq --out " + outputPriorFilename3IBDStitch + endCommand;

	out("Creating frequency data file to calculate priors using PLINK command:\n"); out(plinkCommand4); out("\n\n");

	plinkErrorCode = system(plinkCommand4.c_str());

	if(plinkErrorCode != 0)
	{
		outErr("Problem executing PLINK (error code "); outErr(plinkErrorCode); outErr(") when creating data file with:\n\n ");
		outErr(plinkCommand4); outErr("\n\n");
		outErr("Log file: "); outErr(outputPriorFilename3IBDStitch); outErr(".log:\n\n");
		string plinkLog = "more " + outputPriorFilename3IBDStitch + ".log";
		system(plinkLog.c_str());

		deleteAllTemporaryFiles();
		exit(plinkErrorCode);
	};

	//create the .markers file for ibd_stitch
	string markersFile = outputPriorFilename3IBDStitch + ".markers";

	string markersCommand = "R --vanilla --args " + outputPriorFilename3IBDStitch + ".ped " + outputPriorFilename3IBDStitch + ".map "
		+ outputPriorFilename3IBDStitch + ".frq " + markersFile + " < createIBDStitchDataLIN.R" + endCommand;

	out("Creating markers data file for ibd_stitch with command:\n"); out(markersCommand); out("\n\n");

	int rErrorCode = system(markersCommand.c_str());

	if(rErrorCode != 0)
	{
		outErr("Problem executing R (error code "); outErr(rErrorCode); outErr(") when creating markers data file with:\n\n ");
		outErr(markersCommand); outErr("\n\n");
		
		deleteAllTemporaryFiles();
		exit(rErrorCode);
	};

	//get .par file for ibd_stitch
	string paraFilename = "tempRIA-paras" + toString(pid) + ".par";
	string parCommand = "cp parasTemplate.par " + paraFilename + endCommand;

	out("Copying command:\n"); out(parCommand); out("\n\n");

	unsigned int parErrorCode = system(parCommand.c_str());

	if(parErrorCode != 0)
	{
		outErr("Problem executing copy (error code "); outErr(parErrorCode); outErr(") when copying with:\n\n ");
		outErr(parCommand); outErr("\n\n");

		deleteAllTemporaryFiles();
		exit(parErrorCode);
	};

	//replace strings in .par file
	string resultsFile = outputPriorFilename3IBDStitch + ".out";

	string replaceCommand = "sed -i 's/MARKERFILE/" + markersFile + "/g' " + paraFilename + endCommand;
	out("Replacing strings in parameter file:\n"); out(replaceCommand); out("\n\n");

	unsigned int repErrorCode = system(replaceCommand.c_str());

	if(repErrorCode != 0)
	{
		outErr("Problem executing sed (error code "); outErr(repErrorCode); outErr(") with:\n\n ");
		outErr(replaceCommand); outErr("\n\n");

		deleteAllTemporaryFiles();
		exit(repErrorCode);
	};

	replaceCommand = "sed -i 's/RESULTSFILE/" + resultsFile + "/g' " + paraFilename + endCommand;
	out("Replacing strings in parameter file:\n"); out(replaceCommand); out("\n\n");
	repErrorCode = system(replaceCommand.c_str());

	if(repErrorCode != 0)
	{
		outErr("Problem executing sed (error code "); outErr(repErrorCode); outErr(") with:\n\n ");
		outErr(replaceCommand); outErr("\n\n");

		deleteAllTemporaryFiles();
		exit(repErrorCode);
	};


	//run ibd_stitch!
	string ibdStitchCommand = ibdStitch + " " + paraFilename + endCommand;
	out("Calculating IBD graphs with ibd_stitch:\n"); out(ibdStitchCommand); out("\n\n");
	unsigned int ibdErrorCode = system(ibdStitchCommand.c_str());

	if(ibdErrorCode != 0)
	{
		outErr("Problem executing ibd_stitch (error code "); outErr(ibdErrorCode); outErr(") with:\n\n ");
		outErr(ibdStitchCommand); outErr("\n\n");

		deleteAllTemporaryFiles();
		exit(repErrorCode);
	};

	//calculate IBD probs from ibd graphs
	string ibdCommand = "R --vanilla --args " + resultsFile + " " + outputPriorFilename2 + ".bim " + outputPriorFilename2 + ".fam " + outputPriorFilename2 + "-ibdSt.dat < calcIBDsPedLin.R" + endCommand;
	out("Estimating IBDs from IBD graphs in R with:\n"); out(ibdCommand); out("\n\n");
	ibdErrorCode = system(ibdStitchCommand.c_str());

	if(ibdErrorCode != 0)
	{
		outErr("Problem calculating IBDs in R (error code "); outErr(ibdErrorCode); outErr(") with:\n\n ");
		outErr(ibdCommand); outErr("\n\n");

		deleteAllTemporaryFiles();
		exit(repErrorCode);
	};
	
	//make copy of .fam file for use with posterior calculations
	string copyFamCommand = cpStr + outputPriorFilename2 + ".fam " + postFamFile + endCommand;
	unsigned int cpErrorCode = system(copyFamCommand.c_str());
	if(cpErrorCode != 0)
	{
		outErr("Problem copying (error code "); outErr(cpErrorCode); outErr(") .fam file for posterior calculations with:\n\n ");
		outErr(copyFamCommand); outErr("\n\n");

		deleteAllTemporaryFiles();
		exit(cpErrorCode);
	};

	//read in data and store priors
	priorFilename = outputPriorFilename2 + "-ibdSt.dat";
	setupPriorsUsingInputFile(); //using above file, which was just created

	//copy prior file if req'd
	if(priorOutputFilename != "")
	{
		//make copy of prior ibds
		string copyCommand = cpStr + priorFilename + " " + priorOutputFilename + endCommand;
		unsigned int cpErrorCode = system(copyCommand.c_str());
		if(cpErrorCode != 0)
		{
			outErr("Problem copying (error code "); outErr(cpErrorCode); outErr(") priors file with:\n\n ");
			outErr(copyCommand); outErr("\n\n");

			deleteAllTemporaryFiles();
			exit(cpErrorCode);
		};
	};

	out("Number of affected relative pairs (ARPs) in priors: "); out(priors->snpIBDEstimates.size()); out("\n\n");

	//ostringstream aStringStream;
	//aStringStream << pid;
	//string pidStr = aStringStream.str();

	if(!priorOnly)
	{
		//just output the posterior command here once for convience, so is not outputted for every analysis SNP
		ibdStitchCommand = ibdStitch + " " + paraFilename + endCommand;

		out("Calculating posteriors (for each SNP window) using ibd_stitch command:\n"); out(ibdStitchCommand); out("\n\n");
	};

	deletePostTemporaryFiles();
};

//! Calculates the posteriors using IBD calculations.
void Analysis::calculatePosteriors(unsigned int & windowNumber)
{
	if(doTruffle) { calculatePosteriorsTruffle(); return; }
	else if(doIbdStitch) { calculatePosteriorsIBDStitch(); return; };

	if(posteriorInputFilePrefix != "") outputPostFilename = posteriorInputFilePrefix + "-" + toStringUI(windowNumber); //use previously calculated posteriors
	else
	{
		if(posteriorOutputFilePrefix != "") outputPostFilename = posteriorOutputFilePrefix + "-" + toStringUI(windowNumber);

		snpWindow->createFilesForPosteriorCalc(outputPostFilename);

		//calculate the IBDs for the pruned data set in order to calculate the priors
		string kingCommand = king + " -b " + outputPostFilename + ".bed --homog --prefix " + outputPostFilename + endCommand;

		unsigned int kingErrorCode = system(kingCommand.c_str());

		if(kingErrorCode != 0)
		{
			outErr("Problem executing KING (error code "); outErr(kingErrorCode); outErr(") when estimating IBDs with:\n\n ");
			outErr(kingCommand); outErr("\n\n");

			deleteAllTemporaryFiles();
			exit(kingErrorCode);
		};
	};

	//read in data and store posteriors
	setupPosteriors(outputPostFilename);
	deletePostTemporaryFiles();
};

//! Calculates the posteriors using IBD calculations given by TRUFFLE.
void Analysis::calculatePosteriorsTruffle()
{
	snpWindow->createFilesForPosteriorCalc(outputPostFilename);

	//firstly convert to a VCF gzipped file for use with truffle
	string plinkCommand3 = plink + " --noweb --bfile " + outputPostFilename + " --recode vcf --out " + outputPostFilename + endCommand;

	unsigned int plinkErrorCode = system(plinkCommand3.c_str());

	if(plinkErrorCode != 0)
	{
		outErr("Problem executing PLINK (error code "); outErr(plinkErrorCode); outErr(") when creating data file with:\n\n ");
		outErr(plinkCommand3); outErr("\n\n");
		outErr("Log file: "); outErr(outputPostFilename); outErr(".log:\n\n");
		string plinkLog = "more " + outputPostFilename + ".log";
		system(plinkLog.c_str());

		deleteAllTemporaryFiles();
		exit(plinkErrorCode);
	};

	//now gzip the vcf file
	string gzipCommand = gzip + " -f " + outputPostFilename + ".vcf" + endCommand;

	unsigned int gzipErrorCode = system(gzipCommand.c_str());

	if (gzipErrorCode != 0)
	{
		outErr("Problem executing gzip (error code "); outErr(gzipErrorCode); outErr(") when gzipping data file with:\n\n ");
		outErr(gzipCommand); outErr("\n\n");

		deleteAllTemporaryFiles();
		exit(gzipErrorCode);
	};

	//calculate the IBDs for the posterior
	string truffleCommand = truffle + " --pairs-file " + trufflePairsFilename + " --vcf " + outputPostFilename + ".vcf.gz " + truffleOptions + " --out " + outputPostFilename + " " + endCommand;

	unsigned int truffleErrorCode = system(truffleCommand.c_str());

	if(truffleErrorCode != 0)
	{
		outErr("Problem executing TRUFFLE (error code "); outErr(truffleErrorCode); outErr(") when estimating IBDs with:\n\n ");
		outErr(truffleCommand); outErr("\n\n");

		deleteAllTemporaryFiles();
		exit(truffleErrorCode);
	};

	//read in data and store posteriors
	setupPosteriorsTruffle(outputPostFilename);
	deletePostTemporaryFiles();
};

//! Calculates the posteriors using IBD calculations given by ibd_stitch.
void Analysis::calculatePosteriorsIBDStitch()
{
	snpWindow->createFilesForPosteriorCalc(outputPostFilename);

	//create additioanl files req'd by ibd_stitch
	string prefix = filename.substr(0, filename.length() - 4); //original data file

	//make .ped file this alleles recoded as 1 and 2
	string plinkCommand3 = plink + " --noweb --bfile " + outputPostFilename + " --recode12 --out " + outputPostFilename + endCommand;

	int plinkErrorCode = system(plinkCommand3.c_str());

	if(plinkErrorCode != 0)
	{
		outErr("Problem executing PLINK (error code "); outErr(plinkErrorCode); outErr(") when creating data file with:\n\n ");
		outErr(plinkCommand3); outErr("\n\n");
		outErr("Log file: "); outErr(outputPriorFilename3IBDStitch); outErr(".log:\n\n");
		string plinkLog = "more " + outputPriorFilename3IBDStitch + ".log";
		system(plinkLog.c_str());

		deleteAllTemporaryFiles();
		exit(plinkErrorCode);
	};

	//make .frq file with plink
	string plinkCommand4 = plink + " --noweb --file " + outputPostFilename + " --freq --out " + outputPostFilename + endCommand;

	plinkErrorCode = system(plinkCommand4.c_str());

	if(plinkErrorCode != 0)
	{
		outErr("Problem executing PLINK (error code "); outErr(plinkErrorCode); outErr(") when creating data file with:\n\n ");
		outErr(plinkCommand4); outErr("\n\n");
		outErr("Log file: "); outErr(outputPriorFilename3IBDStitch); outErr(".log:\n\n");
		string plinkLog = "more " + outputPriorFilename3IBDStitch + ".log";
		system(plinkLog.c_str());

		deleteAllTemporaryFiles();
		exit(plinkErrorCode);
	};

	//create the .markers file for ibd_stitch
	string markersFile = outputPostFilename + ".markers";

	string markersCommand = "R --vanilla --args " + outputPostFilename + ".ped " + outputPostFilename + ".map "
		+ outputPostFilename + ".frq " + prefix + ".bim " + markersFile + " < createIBDStitchDataLIN.R" + endCommand;

	int rErrorCode = system(markersCommand.c_str());

	if(rErrorCode != 0)
	{
		outErr("Problem executing R (error code "); outErr(rErrorCode); outErr(") when creating markers data file with:\n\n ");
		outErr(markersCommand); outErr("\n\n");

		deleteAllTemporaryFiles();
		exit(rErrorCode);
	};

	//get .par file for ibd_stitch
	string paraFilename = "tempRIA-paras" + toString(pid) + ".par";
	string parCommand = "cp parasTemplate.par " + paraFilename + endCommand;

	unsigned int parErrorCode = system(parCommand.c_str());

	if(parErrorCode != 0)
	{
		outErr("Problem executing copy (error code "); outErr(parErrorCode); outErr(") when copying with:\n\n ");
		outErr(parCommand); outErr("\n\n");

		deleteAllTemporaryFiles();
		exit(parErrorCode);
	};

	//replace strings in .par file
	string resultsFile = outputPostFilename + ".out";

	string replaceCommand = "sed -i 's/MARKERFILE/" + markersFile + "/g' " + paraFilename + endCommand;

	unsigned int repErrorCode = system(replaceCommand.c_str());

	if(repErrorCode != 0)
	{
		outErr("Problem executing sed (error code "); outErr(repErrorCode); outErr(") with:\n\n ");
		outErr(replaceCommand); outErr("\n\n");

		deleteAllTemporaryFiles();
		exit(repErrorCode);
	};

	replaceCommand = "sed -i 's/RESULTSFILE/" + resultsFile + "/g' " + paraFilename + endCommand;

	repErrorCode = system(replaceCommand.c_str());

	if(repErrorCode != 0)
	{
		outErr("Problem executing sed (error code "); outErr(repErrorCode); outErr(") with:\n\n ");
		outErr(replaceCommand); outErr("\n\n");

		deleteAllTemporaryFiles();
		exit(repErrorCode);
	};

	//run ibd_stitch!
	string ibdStitchCommand = ibdStitch + " " + paraFilename + endCommand;

	unsigned int ibdErrorCode = system(ibdStitchCommand.c_str());

	if(ibdErrorCode != 0)
	{
		outErr("Problem executing ibd_stitch (error code "); outErr(ibdErrorCode); outErr(") with:\n\n ");
		outErr(ibdStitchCommand); outErr("\n\n");

		deleteAllTemporaryFiles();
		exit(repErrorCode);
	};


	//calculate IBD probs from ibd graphs
	string ibdCommand = "R --vanilla --args " + resultsFile + " " + outputPostFilename + ".bim " + outputPostFilename + ".fam " + outputPostFilename + "-ibdSt.dat < calcIBDsPedLin.R" + endCommand;

	ibdErrorCode = system(ibdStitchCommand.c_str());

	if(ibdErrorCode != 0)
	{
		outErr("Problem calculating IBDs in R (error code "); outErr(ibdErrorCode); outErr(") with:\n\n ");
		outErr(ibdCommand); outErr("\n\n");

		deleteAllTemporaryFiles();
		exit(repErrorCode);
	};

	
	//read in data and store posteriors
	ifstream readPostFile(resultsFile.c_str());

	string famID, indivID1, indivID2, arpID;
	bool errorCode;
	double f1j, f2j;

	//read in header
	readPostFile >> famID >> indivID1 >> indivID2 >> famID >> indivID1 >> indivID2;

	//loop thro' priors and store data
	do{
		readPostFile >> famID >> indivID1 >> indivID2;
		if(readPostFile.eof()) break;

		//create string to represent the affected relative pair
		arpID = famID + "-" + indivID1 + "-" + indivID2;

		readPostFile >> f1j >> f2j >> errorCode;

		//priors->snpIBDEstimates[arpID] = new IBDData(f1j, f2j, errorCode);
		posteriors->snpIBDEstimates[arpID] = new IBDData(f1j, f2j, errorCode);

	} while(!readPostFile.eof());

	readPostFile.close();

	deletePostTemporaryFiles();;
};

//! Checks input from KING and deals with nans.
void checkInput(ifstream & readKINGIBDs, bool & hasError, string & badNumber)
{
	if(readKINGIBDs.fail())
	{
		readKINGIBDs.clear();
		readKINGIBDs >> badNumber;
		hasError = true;
	};
};

//! Reads in kinship and IBD estimates and sets up prior data.
void Analysis::setupPriors(string & filenameIBD)
{
	checkKingKinFile(filenameIBD);

	string filenameKINGIBDs = filenameIBD + ".kin";

	ifstream readKINGIBDs;
	readKINGIBDs.open(filenameKINGIBDs.c_str());
	if(!readKINGIBDs.is_open())
	{
		outErr("\nCannot read KING IBD priors: "); outErr(filenameKINGIBDs); outErr("!\n");

		deleteAllTemporaryFiles();
		exit(1);
	};

	//output prior estimates if req'd
	bool outputPriorEst = (priorOutputFilename != "");
	ofstream priorOutFile;
	if(outputPriorEst)
	{
		priorOutFile.open(priorOutputFilename.c_str());
		priorOutFile.precision(10);
		//output header
		priorOutFile << "FAMID INDIV1 INDIV2 IBD1 IBD2 ERROR\n";
	};

	string famID, indivID1, indivID2;
	string prevFamID = "", prevIndivID1 = "", prevIndivID2 = "";
	unsigned int noSNPs;
	double z0, phi, IBD0 = 0, kinship = 0, errorCode;
	double sumf1jf2j, f1j, f2j, sumAll;
	string badNumber, dummy;
	bool hasError;
	string arpID;

	bool errorCodeIsTrue;

	//read in header
	//readKINGIBDs >> famID >> indivID1 >> indivID2 >> famID >> indivID1 >> indivID2 >> famID >> indivID1 >> indivID2;

	for(unsigned int col = 1; col <= noColumns; ++col) readKINGIBDs >> famID;

	//loop thro' subjects and store data
	do{		
		hasError = false;

		//read stuff in
		for(unsigned int col = 1; col <= noColumns; ++col)
		{
			if(col == famIDCol) readKINGIBDs >> famID;
			else if(col == indivID1Col) readKINGIBDs >> indivID1;
			else if(col == indivID2Col) readKINGIBDs >> indivID2;
			else if(col == noSNPsCol) readKINGIBDs >> noSNPs;
			else if(col == z0Col) {readKINGIBDs >> z0; checkInput(readKINGIBDs, hasError, badNumber);}
			else if(col == phiCol) {readKINGIBDs >> phi; checkInput(readKINGIBDs, hasError, badNumber);}
			else if(col == IBD0Col) {readKINGIBDs >> IBD0; checkInput(readKINGIBDs, hasError, badNumber);}
			else if(col == kinshipCol) {readKINGIBDs >> kinship; checkInput(readKINGIBDs, hasError, badNumber);}
			else if(col == errorCodeCol) {readKINGIBDs >> errorCode; checkInput(readKINGIBDs, hasError, badNumber);}
			else readKINGIBDs >> dummy;
		}; 

		errorCodeIsTrue = false; //(errorCode != 0); will be different to expected if calc'd on interval with effect

		if(!(prevFamID == famID && prevIndivID1 == indivID1 && prevIndivID2 == indivID2))
		{
			//create string to represent the affected relative pair
			arpID = famID + "-" + indivID1 + "-" + indivID2;

			//check if data is ok, and if so convert to IBD estimates	
			if(hasError)
			{
				priors->snpIBDEstimates[arpID] = new IBDData(0, 0, true); //record error, bad data
			}
			else
			{
				//store estimate data
				//constrain kinship and IBD0 firstly
				if(kinship < 0) kinship = 0; else if(kinship > 0.5) kinship = 0.5;
				if(IBD0 < 0) IBD0 = 0; else if(IBD0 > 1) IBD0 = 1;

				sumf1jf2j = 1 - IBD0;

				f2j = kinship*4 - sumf1jf2j;
				if(f2j < 0) f2j = 0; else if(f2j > 1) f2j = 1;

				f1j = sumf1jf2j - f2j;
				if(f1j < 0) f1j = 0; else if(f1j > 1) f1j = 1;
			
				sumAll = IBD0 + f1j + f2j;
				f1j = f1j/sumAll;
				f2j = f2j/sumAll;

				//results estimates		
				priors->snpIBDEstimates[arpID] = new IBDData(f1j, f2j, errorCodeIsTrue);
			};
		
			//output prior estimates if req'd
			if(outputPriorEst)
			{
				if(hasError) priorOutFile << famID << " " << indivID1 << " " << indivID2 << " 0 0 1\n";
				else priorOutFile << famID << " " << indivID1 << " " << indivID2 << " " << f1j << " " << f2j << " " << errorCodeIsTrue << "\n";
			};

		};

		prevFamID = famID; prevIndivID1 = indivID1; prevIndivID2 = indivID2;
	}while(!readKINGIBDs.eof());

	readKINGIBDs.close();

	if(outputPriorEst) priorOutFile.close();		
};

//! Separates family and individual IDs.
void getFamAndIndivIDs(string & ID1, string & ID2, string & famID1, string & famID2, string & indivID1, string & indivID2)
{
	unsigned int pos = ID1.find("_");
	famID1 = ID1.substr(0, pos);
	indivID1 = ID1.substr((pos+1));

	pos = ID2.find("_");
	famID2 = ID2.substr(0, pos);
	indivID2 = ID2.substr((pos + 1));

	unsigned int pos1 = indivID1.find("_");
	unsigned int pos2 = indivID2.find("_");
	
	if(pos1 < indivID1.length() || pos2 < indivID2.length())
	{
		outErr("Do not use underscores for IDs when using TRUFFLE (as it uses VCF files)!\n");
		outErr(ID1); outErr(""); outErr(ID2); outErr("\n");
		exit(1);
	};
	
};

//! Reads in kinship and IBD estimates and sets up prior data from TRUFFLE ouput.
void Analysis::setupPriorsTruffle(string & filenameIBD)
{	
	string filenameTruffleIBDs = filenameIBD + ".ibd";

	ifstream readTruffleIBDs;
	readTruffleIBDs.open(filenameTruffleIBDs.c_str());
	if(!readTruffleIBDs.is_open())
	{
		outErr("\nCannot read TRUFFLE IBD priors: "); outErr(filenameTruffleIBDs); outErr("!\n");

		deleteAllTemporaryFiles();
		exit(1);
	};

	//output prior estimates if req'd
	bool outputPriorEst = (priorOutputFilename != "");
	ofstream priorOutFile;
	if(outputPriorEst)
	{
		priorOutFile.open(priorOutputFilename.c_str());
		priorOutFile.precision(10);
		//output header
		priorOutFile << "FAMID INDIV1 INDIV2 IBD1 IBD2 ERROR\n";
	};

	string ID1, ID2;
	string famID1, famID2, indivID1, indivID2;
	string prevFamID, prevIndivID1, prevIndivID2;
	double IBD0 = 0, IBD1 = 0, IBD2 = 0;
	string nMark, nCommon, ibd1Max, ibd1nSegs, ibd2Max, ibd2nSegs, sex;
	noColumns = 12;

	bool hasError;
	string arpID;

	bool errorCodeIsTrue = false;

	//read in header
	//    ID1   ID2    NMARK  NCOMMON         IBD0        IBD1_MAX IBD1_NSEGS      IBD1        IBD2_MAX IBD2_NSEGS      IBD2   SEX
	//

	for (unsigned int col = 1; col <= noColumns; ++col) readTruffleIBDs >> ID1;

	//loop thro' subjects and store data
	do{
		hasError = false;

		//read stuff in		
		readTruffleIBDs >> ID1 >> ID2 >> nMark >> nCommon >> IBD0 >> ibd1Max >> ibd1nSegs >> IBD1 >> ibd2Max >> ibd2nSegs >> IBD2 >> sex;		
		
		getFamAndIndivIDs(ID1, ID2, famID1, famID2, indivID1, indivID2);

		if(famID1 == famID2 && !(prevFamID == famID1 && prevIndivID1 == indivID1 && prevIndivID2 == indivID2))
		{
			//create string to represent the affected relative pair
			arpID = famID1 + "-" + indivID1 + "-" + indivID2;

			//check if data is ok, and if so convert to IBD estimates	
			if(!hasError)
			{
				//results estimates		
				priors->snpIBDEstimates[arpID] = new IBDData(IBD1, IBD2, errorCodeIsTrue);
			}
			else
			{				
				priors->snpIBDEstimates[arpID] = new IBDData(0, 0, true); //record error, bad data
			};

			//output prior estimates if req'd
			if (outputPriorEst)
			{
				if(!hasError) priorOutFile << famID1 << " " << indivID1 << " " << indivID2 << " " << IBD1 << " " << IBD2 << " " << errorCodeIsTrue << "\n"; 
				else priorOutFile << famID1 << " " << indivID1 << " " << indivID2 << " 0 0 1\n";
			};

		};

		prevFamID = famID1; prevIndivID1 = indivID1; prevIndivID2 = indivID2;
	}while(!readTruffleIBDs.eof());

	readTruffleIBDs.close();

	if (outputPriorEst) priorOutFile.close();
};

//! Outputs prior IBD estimates based on whole SNP data.
void Analysis::setupPriorsUsingInputFile()
{
	ifstream readPriorFile(priorFilename.c_str());

	if(!readPriorFile.is_open())
	{
		outErr("\nCannot read prior file: "); outErr(priorFilename); outErr("!\n");
		deleteAllTemporaryFiles();
		exit(1);
	};

	string famID, indivID1, indivID2, arpID;
	bool errorCode;
	double f1j, f2j;

	//read in header
	readPriorFile >> famID >> indivID1 >> indivID2 >> famID >> indivID1 >> indivID2;

	//loop thro' priors and store data
	do{
		readPriorFile >> famID >> indivID1 >> indivID2;
		if(readPriorFile.eof()) break;

		//create string to represent the affected relative pair
		arpID = famID + "-" + indivID1 + "-" + indivID2;

		readPriorFile >> f1j >> f2j >> errorCode;

		priors->snpIBDEstimates[arpID] = new IBDData(f1j, f2j, errorCode);

	}while(!readPriorFile.eof());

	readPriorFile.close();

	ostringstream aStringStream;
	aStringStream << pid;
	string pidStr = aStringStream.str();

	out("Number of affected relative pairs (ARPs) in priors: "); out(priors->snpIBDEstimates.size()); out("\n\n");

	//just output the posterior command here once for convience, so is not outputted for every analysis SNP
	if(doTruffle)
	{
		string truffleCommand = truffle + " --vcf " + outputPostFilename + ".vcf.gz " + truffleOptions + " --out" + outputPostFilename + " " + endCommand;
		out("Calculating posteriors (for each SNP window) using TRUFFLE command:\n"); out(truffleCommand); out("\n\n");
	}
	else if(doIbdStitch)
	{
		string ibdStitchCommand = ibdStitch + " " + "paras" + toString(pid) + ".par" + endCommand;
		out("Calculating posteriors (for each SNP window) using ibd_stitch command:\n"); out(ibdStitchCommand); out("\n\n");
	}
	else
	{
		string kingCommand = king + " -b " + outputPostFilename + ".bed --homog --prefix " + outputPostFilename + " " + endCommand;
		out("Calculating posteriors (for each SNP window) using KING command:\n"); out(kingCommand); out("\n\n");
	};
};



//! Checks header of KING for posterior files.
void Analysis::checkKingKinFile(string & filenameIBD)
{
	string filenameKINGIBDs = filenameIBD + ".kin";

	ifstream readKINGIBDs;
	readKINGIBDs.open(filenameKINGIBDs.c_str());
	if(!readKINGIBDs.is_open())
	{
		outErr("\nCannot read KING IBD posteriors: "); outErr(filenameKINGIBDs); outErr("!\n");

		deleteAllTemporaryFiles();
		exit(1);
	};

	char returnChar = '\n';
	set<string> citeLines;
	string aLine;
	string aChar;
	list<string> words;

	getline(readKINGIBDs, aLine, returnChar);

	unsigned int len = aLine.length();
	unsigned int prevPos = 0;

	for(unsigned int i = 0; i < len; ++i)
	{
		aChar = aLine.substr(i, 1);
		if(aChar == " " || aChar == "\t") 
		{
			words.push_back(aLine.substr(prevPos, i - prevPos));
			prevPos = i + 1;			
		}
		else if(i == len-1) 
		{
			words.push_back(aLine.substr(prevPos, i - prevPos + 1));
		};
	};

	readKINGIBDs.close();


	noColumns = words.size();
	//FID     ID1     ID2     N_SNP   Z0      Phi     IBD0    Kinship Error
	famIDCol = 0; indivID1Col = 0; indivID2Col = 0; noSNPsCol = 0; z0Col = 0; phiCol = 0; IBD0Col = 0; kinshipCol = 0; errorCodeCol = 0;

	unsigned int colNo = 1;
	for(list<string>::const_iterator w = words.begin(); w != words.end(); ++w, ++colNo)
	{
		if(*w == "FID") famIDCol = colNo;
		else if(*w == "ID1") indivID1Col = colNo;
		else if(*w == "ID2") indivID2Col = colNo;
		else if(*w == "N_SNP") noSNPsCol = colNo;
		else if(*w == "Z0") z0Col = colNo;
		else if(*w == "Phi") phiCol = colNo;
		else if(*w == "IBD0") IBD0Col = colNo;
		else if(*w == "Kinship") kinshipCol = colNo;
		else if(*w == "Error") errorCodeCol = colNo;
	};

	if(famIDCol == 0 || indivID1Col  == 0 || indivID2Col  == 0 || noSNPsCol  == 0 || z0Col  == 0 || phiCol  == 0 || IBD0Col  == 0 || kinshipCol  == 0 || errorCodeCol  == 0)
	{
		
		outErr("The following columns are required in the KING IBD output!\n");
		outErr("FID, ID1, ID2, N_SNP, Z0, Phi, IBD0, Kinship and Error!\n");
		outErr("You may need to use an older version of KING if the output files have changed.\n");
		deleteAllTemporaryFiles();
		exit(1);
	};
	
};


//! Reads in kinship and IBD estimates and sets up posterior data.
void Analysis::setupPosteriors(string & filenameIBD)
{
	checkKingKinFile(filenameIBD);
	string filenameKINGIBDs = filenameIBD + ".kin";

	ifstream readKINGIBDs;
	readKINGIBDs.open(filenameKINGIBDs.c_str());
	if(!readKINGIBDs.is_open())
	{
		outErr("\nCannot read KING IBD posteriors: "); outErr(filenameKINGIBDs); outErr("!\n");

		deleteAllTemporaryFiles();
		exit(1);
	};
	string famID, indivID1, indivID2;
	string prevFamID = "", prevIndivID1 = "", prevIndivID2 = "";
	unsigned int noSNPs;
	double z0, phi, IBD0 = 0, kinship = 0, errorCode;
	double sumf1hatjf2hatj, f1hatj, f2hatj, sumAll;
	string badNumber, dummy;
	bool hasError;
	string arpID;	

	bool errorCodeIsTrue;

	//read in header
	//readKINGIBDs >> famID >> indivID1 >> indivID2 >> famID >> indivID1 >> indivID2 >> famID >> indivID1 >> indivID2;
	for(unsigned int col = 1; col <= noColumns; ++col) readKINGIBDs >> famID;
	//loop thro' priors and set corresponding posterior estimates to zero if prior is zero
	map<string, IBDData *>::const_iterator pri = priors->snpIBDEstimates.begin();

	//loop thro' subject pairs and store data
	do{		
		hasError = false;
		
		//read stuff in
		for(unsigned int col = 1; col <= noColumns; ++col)
		{
			if(col == famIDCol) readKINGIBDs >> famID;
			else if(col == indivID1Col) readKINGIBDs >> indivID1;
			else if(col == indivID2Col) readKINGIBDs >> indivID2;
			else if(col == noSNPsCol) readKINGIBDs >> noSNPs;
			else if(col == z0Col) {readKINGIBDs >> z0; checkInput(readKINGIBDs, hasError, badNumber);}
			else if(col == phiCol) {readKINGIBDs >> phi; checkInput(readKINGIBDs, hasError, badNumber);}
			else if(col == IBD0Col) {readKINGIBDs >> IBD0; checkInput(readKINGIBDs, hasError, badNumber);}
			else if(col == kinshipCol) {readKINGIBDs >> kinship; checkInput(readKINGIBDs, hasError, badNumber);}
			else if(col == errorCodeCol) {readKINGIBDs >> errorCode; checkInput(readKINGIBDs, hasError, badNumber);}
			else readKINGIBDs >> dummy;
		}; 

		errorCodeIsTrue = false; //(errorCode != 0); will be different to expected if calc'd on interval with effect

		if(!(prevFamID == famID && prevIndivID1 == indivID1 && prevIndivID2 == indivID2))
		{

			//create string to represent the affected relative pair
			arpID = famID + "-" + indivID1 + "-" + indivID2;

			//do not duplicate the last row		
			if(hasError)
			{
				posteriors->snpIBDEstimates[arpID] = new IBDData(0, 0, true); //record error, bad data

				if(pri == priors->snpIBDEstimates.end() || pri->first != arpID)
				{
					//find correct pair in posteriors
					pri = priors->snpIBDEstimates.find(arpID);			
				};
			}
			else
			{
				//store estimate data
				if(kinship < 0) kinship = 0; else if(kinship > 0.5) kinship = 0.5;
				if(IBD0 < 0) IBD0 = 0; else if(IBD0 > 1) IBD0 = 1;

				sumf1hatjf2hatj = 1 - IBD0;
				f2hatj = kinship*4 - sumf1hatjf2hatj;
				if(f2hatj < 0) f2hatj = 0; else if(f2hatj > 1) f2hatj = 1;

				f1hatj = sumf1hatjf2hatj - f2hatj;
				if(f1hatj < 0) f1hatj = 0; else if(f1hatj > 1) f1hatj = 1;
	
				//check prior and posteriors are for the same pair
				if(pri == priors->snpIBDEstimates.end() || pri->first != arpID)
				{
					//find correct pair in posteriors
					pri = priors->snpIBDEstimates.find(arpID);			
				};

				//set posteriors to zero if priors are zero
				if(pri != priors->snpIBDEstimates.end())
				{		
					if(pri->second->f1j == 0) f1hatj = 0;
					if(pri->second->f2j == 0) f2hatj = 0;
					if((pri->second->f1j + pri->second->f2j) == 1) IBD0 = 0;
				};

				//rescale to sum to 1
				sumAll = IBD0 + f1hatj + f2hatj;

				//could be all zero if set to zero from prior, make sure there is no problem 
				if(sumAll > 0)
				{
					f1hatj = f1hatj/sumAll;
					f2hatj = f2hatj/sumAll;

					//results estimates		
					posteriors->snpIBDEstimates[arpID] = new IBDData(f1hatj, f2hatj, errorCodeIsTrue);
				}
				else
				{
					posteriors->snpIBDEstimates[arpID] = new IBDData(f1hatj, f2hatj, true); //error
				};
			};
		};
		
		prevFamID = famID; prevIndivID1 = indivID1; prevIndivID2 = indivID2;
		++pri;
	}while(!readKINGIBDs.eof());
	
	readKINGIBDs.close();
};

//! Reads in kinship and IBD estimates and sets up posterior data using Truffle.
void Analysis::setupPosteriorsTruffle(string & filenameIBD)
{
	string filenameTruffleIBDs = filenameIBD + ".ibd";

	ifstream readTruffleIBDs;
	readTruffleIBDs.open(filenameTruffleIBDs.c_str());
	if (!readTruffleIBDs.is_open())
	{
		outErr("\nCannot read TRUFFLE IBD posteriors: "); outErr(filenameTruffleIBDs); outErr("!\n");

		deleteAllTemporaryFiles();
		exit(1);
	};

	string ID1, ID2;
	string famID1, famID2, indivID1, indivID2;
	string prevFamID, prevIndivID1, prevIndivID2;
	double IBD0 = 0, IBD1 = 0, IBD2 = 0;
	string nMark, nCommon, ibd1Max, ibd1nSegs, ibd2Max, ibd2nSegs, sex;
	noColumns = 12;

	string arpID;

	bool errorCodeIsTrue = false;

	//read in header
	//    ID1   ID2    NMARK  NCOMMON         IBD0        IBD1_MAX IBD1_NSEGS      IBD1        IBD2_MAX IBD2_NSEGS      IBD2   SEX
	//

	for (unsigned int col = 1; col <= noColumns; ++col) readTruffleIBDs >> ID1;

	//loop thro' priors and set corresponding posterior estimates to zero if prior is zero
	map<string, IBDData*>::const_iterator pri;// = priors->snpIBDEstimates.begin();

	//loop thro' subject pairs and store data
	do {		
		//read stuff in
		readTruffleIBDs >> ID1 >> ID2 >> nMark >> nCommon >> IBD0 >> ibd1Max >> ibd1nSegs >> IBD1 >> ibd2Max >> ibd2nSegs >> IBD2 >> sex;

		getFamAndIndivIDs(ID1, ID2, famID1, famID2, indivID1, indivID2);

		//do not duplicate the last row	
		if(famID1 == famID2 && !(prevFamID == famID1 && prevIndivID1 == indivID1 && prevIndivID2 == indivID2))
		{

			//create string to represent the affected relative pair
			arpID = famID1 + "-" + indivID1 + "-" + indivID2;

			pri = priors->snpIBDEstimates.find(arpID);	
			
			if(pri == priors->snpIBDEstimates.end() || pri->second->error)
			{
				posteriors->snpIBDEstimates[arpID] = new IBDData(0, 0, true); //record error, bad data
			}
			else
			{				
				//result estimates		
				posteriors->snpIBDEstimates[arpID] = new IBDData(IBD1, IBD2, errorCodeIsTrue);				
			};
		};

		prevFamID = famID1; prevIndivID1 = indivID1; prevIndivID2 = indivID2;
		
	}while (!readTruffleIBDs.eof());

	readTruffleIBDs.close();
};

//! Reads in SNP data and converts to the number of minor alleles.
unsigned int Analysis::getNextNoOfMinorAlleles(ifstream & readSNPData, unsigned int & bitCount)
{
	int allele1, allele2;
	unsigned int noMinorAlleles = 0;
	int one = '\1';
		
	//read in the next piece of data
	if(bitCount == 9)
	{
		
		readSNPData.read(oneBuffer, 1);
		if(readSNPData.eof())
		{			
			outErr("Error: binary SNP file (.bed) is incomplete!\n");
			deleteAllTemporaryFiles();
			exit(1);
		};
		
		aBit = oneBuffer[0];
			
		bitCount = 1;
	};

	allele1 = aBit & one; //read the least significant bit				
	aBit = aBit >> 1; //shift bits to the right
	allele2 = aBit & one; //read the new least significant bit				
	aBit = aBit >> 1; //shift bits to the right for next time

	bitCount += 2;	

	//if genotype is encoded 1/0 then the genotype is missing so do not add it
	if(allele1 == 1 && allele2 == 1)
	{	
		noMinorAlleles = 0;
	}
	else if(allele1 == 0 && allele2 == 1)
	{	
		noMinorAlleles = 1;
	}
	else if(allele1 == 0 && allele2 == 0)
	{	
		noMinorAlleles = 2;
	}
	else
		noMinorAlleles = 3; //denotes missing genotype

	return noMinorAlleles;
};

//! Calculates the p-value from a Chi square value with 1 df.
double Analysis::getPvalueChiSq1DF(double & chisq)
{	
	double a = sqrt(chisq)*oneOverSqRoot2;
	int ind = 0;
	return erfc1(&ind, &a);
};

//! Fits the model parameters, varA and varD here.
void Analysis::fitModels(const unsigned int & snpID)
{
	FindFit findFit(modelIBD); 	

	//get alt model neg log like, base 10 ready for LOD score
	//no need to calculate null like as constants for data prob cancel out when calc'ing the MLS, maximum likelihood LOD score
	//which is s.t. MLS = -negLog10Alt
	double negLog10Alt;
	double MLS, eval2;
	bool fittedOK;
	bool fittedOK2;
	double paraVarA, paraVarD;
	unsigned int paraToFit = 1;
	
	if(noDomVar)
	{
		unsigned int paraFit = 1;
		modelIBD->setParameter(1, prevVarA);
		modelIBD->setParameter(2, 0);
		
		fittedOK = findFit.newtonsMethodForOneVariable(negLog10Alt, paraFit);
		
		//try fitting with 0
		if(!fittedOK)
		{
			modelIBD->setParameter(1, 0);
			fittedOK = findFit.newtonsMethodForOneVariable(negLog10Alt, paraFit);
		};
		
		if(fittedOK)
		{
			prevVarA = modelIBD->getParameter(1);
		};

		if(modelIBD->getParameter(1) < 0) //variance of varA cannot be <0
		{
			//try large starting value to see if a lower minimum somewhere
			modelIBD->setParameter(1, 10);
			fittedOK2 = findFit.newtonsMethodForOneVariable(eval2, paraFit);

			if(fittedOK2 && modelIBD->getParameter(1) > 0 && eval2 < negLog10Alt)
			{
				MLS = -eval2;
			}
			else
			{
				modelIBD->setParameter(1, 0);
				MLS = 0;
			};
		}
		else
		{
			MLS = -negLog10Alt;
		};

	}
	else
	{

		map<unsigned int, double> parameters;
		set<unsigned int> parasToFit;
		parasToFit.insert(1);
		parasToFit.insert(2);

		//initial default parameters, set VarA and VarD to 0
		modelIBD->setParameter(1, prevVarA);
		modelIBD->setParameter(2, prevVarD);

		//fit the model	
		fittedOK = findFit.newtonsMethod(negLog10Alt, parasToFit);

		//try fitting with downhill simplex if failed to fit or negative parameters, then newtons
		if(!fittedOK || modelIBD->getParameter(1) < 0 || modelIBD->getParameter(2) < 0)
		{
			modelIBD->setParameter(1, 0);
			modelIBD->setParameter(2, 0);
			fittedOK2 = findFit.downhillSimplex(negLog10Alt, parasToFit, 1e-6, 100, 0.01);

			paraVarA = modelIBD->getParameter(1);
			paraVarD = modelIBD->getParameter(2);

			if(fittedOK2) fittedOK = findFit.newtonsMethod(negLog10Alt, parasToFit);

			//if still not fitted, fit accuarily using downhill simplex only
			if(fittedOK2 && (!fittedOK || modelIBD->getParameter(1) < 0 || modelIBD->getParameter(2) < 0))
			{
				modelIBD->setParameter(1, paraVarA);
				modelIBD->setParameter(2, paraVarD);
				fittedOK = findFit.downhillSimplex(negLog10Alt, parasToFit, 1e-8, 1000, 1e-6);

				paraVarA = modelIBD->getParameter(1);
				paraVarD = modelIBD->getParameter(2);

				//try and get a bit more accurate if one parameter is 0
				if(fittedOK && paraVarD == 0)
				{
					paraToFit = 1;
					modelIBD->setParameter(1, paraVarA);
					modelIBD->setParameter(2, 0);					
					fittedOK2 = findFit.newtonsMethodForOneVariable(eval2, paraToFit);

					if(fittedOK2 && modelIBD->getParameter(1) < 0)
					{
						modelIBD->setParameter(1, 0);
						modelIBD->setParameter(2, 0);
						MLS = 0;
					}
					else if(fittedOK2)
					{
						negLog10Alt = eval2;
					}
					else
					{
						modelIBD->setParameter(1, paraVarA); //switch parameters back
						modelIBD->setParameter(2, paraVarD);
					};
				}
				else if(fittedOK && paraVarA == 0)
				{
					paraToFit = 2;
					modelIBD->setParameter(1, 0);
					modelIBD->setParameter(2, paraVarD);					
					fittedOK2 = findFit.newtonsMethodForOneVariable(eval2, paraToFit);

					if(fittedOK2 && modelIBD->getParameter(2) < 0)
					{
						modelIBD->setParameter(1, 0);
						modelIBD->setParameter(2, 0);
						MLS = 0;
					}
					else if(fittedOK2)
					{
						negLog10Alt = eval2;
					}
					else
					{
						modelIBD->setParameter(1, paraVarA); //switch parameters back
						modelIBD->setParameter(2, paraVarD);
					};
				};
			};//end of accurate simplex

		};//end of simplex + newton's

		//set parameters for starting point for next time
		if(fittedOK)
		{
			MLS = -negLog10Alt; 
			//in case converges on answer worse than null parameters
			if(MLS < 0)
			{
				modelIBD->setParameter(1, 0);
				modelIBD->setParameter(2, 0);
				MLS = 0;
			};

			if(MLS==0) MLS = 0; //this looks pretty pointless? but stops the annoying output "-0"
			prevVarA = modelIBD->getParameter(1);
			prevVarD = modelIBD->getParameter(2);

		};
	
	};//end of fitting both vars

	resultsFile << snpID << " " << snpWindow->getAnalysisChromosome() << " " << snpWindow->getAnalysisSNPName()  << " " << snpWindow->getAnalysisSNPCM() << " " << snpWindow->getAnalysisSNPBP() << " "; 
	
	if(fittedOK)
	{
		resultsFile << modelIBD->getParameter(1) << " " << modelIBD->getParameter(2) << " " << MLS << "\n";
	}
	else resultsFile << "NA NA NA\n"; 

};
