/************************************************************************
 * RIA, version 1.00
 * Copyright 2015
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


/*! \file main.cpp
    \brief This file reads in the initial input files and options.
    
    This file also outputs usage instructions and program details.
*/

#include <iostream>
#include <ostream>
#include <string>
#include <time.h>
#include <cstdlib>

using namespace std; // initiates the "std" or "standard" namespace
 
#include "main.h"
#include "Analysis.h"

bool outputToScreen = true; 
ofstream logFile;

//! Output program title to screen
void header()
{
	out("\nRIA: Regional IBD Analysis, v1.00\n");
	out("------------------------------------------------------------\n");
	out("Copyright 2015 Richard Howey, GNU General Public License, v3\n");
	out("Institute of Genetic Medicine, Newcastle University\n\n");
};

//! Outputs program usage to screen.
void usage()
{
		header();
	 	
		out("Usage:\n  ./ria [options] -i pedigree.bed \n");
		out(" or ./ria -pf parameterfile\n\n");

		out("Options:\n");
		out("  -window-size-cm n      -- set window size to n centimorgans\n");
		//out("  -missing-cm-zero       -- treat 0 cM as missing, otherwise only negative values\n");
		out("  -decreasing-cm         -- set decreasing centimorgan values as previous valid value\n");
		out("  -window-min-snps m     -- set required minimium number of SNPs in a window to m\n");
		out("  -step-size s           -- step size of windows, s\n");
		out("  -start-snp a           -- start analysis from SNP number a\n");
		out("  -end-snp b             -- end analysis at SNP number b\n");
		out("  -start-snp-name x      -- start analysis from SNP name x\n");
		out("  -end-snp-name y        -- end analysis at SNP name y\n");
		out("  -job d t               -- job number d of t\n");
		out("  -i file.bed            -- input binary pedigree file, file.bed\n");		
		out("  -o results.dat         -- output results file, results.dat\n");
		out("  -i-prior file          -- input prior IBDs, file\n");
		out("  -o-prior file          -- output prior IBDs, file\n");
		out("  -prior-only            -- calculate prior IBDs only\n");
		out("  -plink command         -- command used to run PLINK\n");
		out("  -plink-options \"ops\"   -- PLINK pruning options used to calculate the prior\n");
		out("  -king command          -- command used to run KING\n");
		out("  -truffle command       -- command used to run Truffle instead of KING\n");
		out("  -truffle-options \"ops\" -- truffle options, e.g. \"--cpu 16\"  \n");
		out("  -gzip command          -- command used to run gzip, needed for use with Truffle\n");
		out("  -log results.log       -- log filename, results.log\n");	
		out("  -ndv                   -- no dominance variance\n");
		out("  -so                    -- suppress output to screen\n\n");

		out("Default Options in Effect:\n");
		out("  -window-size-cm 15\n");
		out("  -window-min-snps 100\n");
		out("  -step-size 50\n");
		out("  -plink plink\n");
		out("  -king king\n");
		out("  -o riaResults.dat\n\n");
};


//! Gets an option value either from the parameter file or the command line.
void getOptionValue(unsigned int & anUnInt, bool & useParaFile, int & argcount, int & argc, char * argv[], ifstream & readParaFile)
{	
	if(useParaFile)
	{
		if(readParaFile.eof()) return;
		readParaFile >> anUnInt;		
	}
	else
	{
		argcount++; if(argcount >= argc) return;				
		anUnInt = atoi(argv[argcount]);		
	};
};

//! Gets an option value either from the parameter file or the command line.
void getOptionValue(double & aDouble, bool & useParaFile, int & argcount, int & argc, char * argv[], ifstream & readParaFile)
{
	if(useParaFile)
	{
		if(readParaFile.eof()) return;		
		readParaFile >> aDouble;		
	}
	else
	{
		argcount++; if(argcount >= argc) return;			
		aDouble = atof(argv[argcount]);		
	};
};

//! Gets an option value either from the parameter file or the command line.
void getOptionValue(string & aString, bool & useParaFile, int & argcount, int & argc, char * argv[], ifstream & readParaFile)
{
	if(useParaFile)
	{
		if(readParaFile.eof()) return;		
		readParaFile >> aString;
	}
	else
	{
		argcount++; if(argcount >= argc) return;		
		aString = argv[argcount];
	};
};

//! Gets the log filename from the results file name.
string getDefaultLogFileName(string & outFileName)
{
	unsigned int filenameLength = outFileName.length();
	string logFileName;

	//find extension
	unsigned int a = filenameLength - 1;
	while(a > 0)
	{
		if(outFileName.substr(a, 1) == ".") break;
		a--;
	};

	if(a > 0) logFileName = outFileName.substr(0, a) + ".log";
	else  logFileName = outFileName + ".log";

	return logFileName;
};

//! The start of the program.
int main(int argc, char * argv[])
{
	time_t start,end;
	double dif;
	time(&start);

	int argcount = 1;
	string option;
	string filename = "";
	string paraFilename = "";
	string outputFilename = "riaResults.dat";
	string logFilename = "";
	string priorFilename = "";
	string priorOutputFilename = "";
	string plink = "plink";
	string king = "king";
	string truffle = "";
	string truffleOptions = "--cpu 4";
	string gzip = "gzip";
	string ibdStitch = "";
	string startSNPName = "";
	string endSNPName = "";
	double windowCMSize = 15;
	bool missingCMZero = false;
	bool decreasingCM = false;
	unsigned int windowStepSize = 50;
	unsigned int windowMinSNPSize = 100;
	unsigned int startSNP = 0;
	unsigned int endSNP = 0;
	unsigned int jobNo = 0;
	unsigned int jobTotal = 0;
	bool priorOnly = false;
	bool noDomVar = false;
	bool keepTempFiles = false;
	string plinkPriorOptions;

	outputToScreen = true;	

	bool useParaFile = false;
	ifstream readParaFile;

	if(argcount < argc) option = argv[argcount];

	//deal with parameter file
	if(option == "-pf")
	{
		argcount++; 
		if(argcount < argc) paraFilename = argv[argcount];

		//open parameter file		
		readParaFile.open(paraFilename.c_str());
		if(!readParaFile.is_open())
		{
			header();
			outErr("Cannot read parameter file: "); outErr(paraFilename); outErr("!\n");			
			exit(1);
		};

		argcount++; 
		useParaFile = true;
	};

	//set given options
	while((!useParaFile && argcount < argc && argv[argcount][0] == '-') || (useParaFile && !readParaFile.eof()))
	{

		if(useParaFile)
		{
			//find the start of the next command
			do{
				readParaFile >> option;
				if(option.length() >= 2 && option.substr(0,1) == "-") break;				
			}while(!readParaFile.eof());
		}
		else
		{
			option = argv[argcount];
		};

		if(useParaFile && readParaFile.eof()) break;

		if(option ==  "-window-size-cm" || option ==  "-window-size-cM" || option ==  "-ws") getOptionValue(windowCMSize, useParaFile, argcount, argc, argv, readParaFile);
		else if(option ==  "-window-min-snps") getOptionValue(windowMinSNPSize, useParaFile, argcount, argc, argv, readParaFile);
		else if(option ==  "-missing-cm") missingCMZero = true;
		else if(option ==  "-decreasing-cm") decreasingCM = true;
		else if(option ==  "-step-size") getOptionValue(windowStepSize, useParaFile, argcount, argc, argv, readParaFile);
		else if(option ==  "-start-snp" || option ==  "-s") getOptionValue(startSNP, useParaFile, argcount, argc, argv, readParaFile);
		else if(option ==  "-end-snp" || option ==  "-e") getOptionValue(endSNP, useParaFile, argcount, argc, argv, readParaFile);
		else if(option ==  "-start-snp-name") getOptionValue(startSNPName, useParaFile, argcount, argc, argv, readParaFile);
		else if(option ==  "-end-snp-name") getOptionValue(endSNPName, useParaFile, argcount, argc, argv, readParaFile);
		else if(option ==  "-i") getOptionValue(filename, useParaFile, argcount, argc, argv, readParaFile);
		else if(option ==  "-i-prior") getOptionValue(priorFilename, useParaFile, argcount, argc, argv, readParaFile);
		else if(option ==  "-o-prior") getOptionValue(priorOutputFilename, useParaFile, argcount, argc, argv, readParaFile);
		else if(option ==  "-plink") getOptionValue(plink, useParaFile, argcount, argc, argv, readParaFile);
		else if(option ==  "-plink-options") getOptionValue(plinkPriorOptions, useParaFile, argcount, argc, argv, readParaFile);
		else if(option ==  "-king") getOptionValue(king, useParaFile, argcount, argc, argv, readParaFile);
		else if(option ==  "-truffle") getOptionValue(truffle, useParaFile, argcount, argc, argv, readParaFile);
		else if(option ==  "-truffle-options") getOptionValue(truffleOptions, useParaFile, argcount, argc, argv, readParaFile);
		else if(option ==  "-gzip") getOptionValue(gzip, useParaFile, argcount, argc, argv, readParaFile);
		else if(option ==  "-log") getOptionValue(logFilename, useParaFile, argcount, argc, argv, readParaFile);	
		else if(option ==  "-ndv") noDomVar = true;
		else if(option ==  "-ibd-stitch") getOptionValue(ibdStitch, useParaFile, argcount, argc, argv, readParaFile);
		else if(option ==  "-keep-temps") keepTempFiles = true;
		else if(option ==  "-prior-only") priorOnly = true;			
		else if(option ==  "-job")
		{	
			getOptionValue(jobNo, useParaFile, argcount, argc, argv, readParaFile);
			getOptionValue(jobTotal, useParaFile, argcount, argc, argv, readParaFile);			
		}
		else if(option ==  "-o") getOptionValue(outputFilename, useParaFile, argcount, argc, argv, readParaFile);
		else if(option == "-so") outputToScreen = false;		
		else
		{
			if(logFilename != "") logFile.open(logFilename.c_str());
			else logFile.open(getDefaultLogFileName(outputFilename).c_str());			

			header();
			if(useParaFile) {outErr("Unrecognised option: "); outErr(option); outErr("\n\n");}
			else {outErr("Unrecognised command line switch: "); outErr(option); outErr("\n\n");};			
    		exit(1);
		};

		if(!useParaFile) argcount++;
	};

	//if(argcount < argc) filename = argv[argcount];	
	
	if(logFilename != "") logFile.open(logFilename.c_str());
	else
	{
		logFilename = getDefaultLogFileName(outputFilename);
		logFile.open(logFilename.c_str());	
	};

	if(filename == "")
	{	
		usage();
		if(argc > 1) outErr("\nInput file not set!\n\n");
		exit(1);
	};

	if (filename.length() >= 4 && filename.substr(filename.length() - 4) != ".bed")
	{
		header();
		outErr("A binary pedigree file (.bed) is required for RIA!\n\n");
		outErr(filename); outErr(" is not a .bed file.\n\n");
		outErr("Try creating a .bed file from a .ped file using PLINK:\n");
		outErr("plink --file mydata --make-bed\n\n");
		exit(1);
	}
	else if (startSNP > endSNP&& endSNP != 0)
	{
		outErr("\nThe start SNP must not be after the end SNP!\n\n");
		exit(1);
	}
	else if (jobNo > jobTotal)
	{
		outErr("\nThe job number cannot be larger than the total number of jobs!\n\n");
		exit(1);
	}
	else if (priorOnly && priorOutputFilename == "")
	{
		outErr("\nYou must set a file name with \"-o-prior\" if you want to output priors to file!\n\n");
		exit(1);
	}
	else if (priorFilename != "" && priorOutputFilename != "")
	{
		outErr("\nYou cannot input and output priors at the same time!\n\n");
		exit(1);
	}
	else if (jobNo != 0 && (startSNP != 0 || endSNP != 0))
	{
		outErr("\nYou cannot set the SNPs to analyse using job numbers and start/end SNPs!\n\n");
		exit(1);
	}
	else if (startSNP != 0 && startSNPName != "")
	{
		outErr("\nYou cannot set the start SNP with a number and a name!\n\n");
		exit(1);
	}
	else if (endSNP != 0 && endSNPName != "")
	{
		outErr("\nYou cannot set the end SNP with a number and a name!\n\n");
		exit(1);
	}
	else if (windowCMSize <= 0)
	{
		outErr("\nThe window size must be strictly positive!\n\n");
		exit(1);
	}
	else if (windowStepSize == 0)
	{
		outErr("\nThe window step size cannot be zero!\n\n");
		exit(1);
	};

	//output options to screen	
	header();
	out("Parameters:\n");
	out("Input file: "); out(filename); out("\n");
	if(!priorOnly) {out("Output file: "); out(outputFilename); out("\n");};

	if(priorOnly) out("Calculating prior IBDs only\n");
	if(priorFilename != "") {out("Prior IBD input file: "); out(priorFilename); out("\n");};
	if(priorOutputFilename != "") {out("Prior IBD output file: "); out(priorOutputFilename); out("\n");};

	out("Log file: "); out(logFilename); out("\n");
	
	if(!priorOnly)
	{
		if(noDomVar) out("Using no dominance variance model\n");
		else {out("Using additive and dominance variance model\n");};
	
		if(startSNP != 0) {out("Start at SNP number: "); out(startSNP); out("\n");}
		else if(startSNPName != "") {out("Start at SNP: "); out(startSNPName); out("\n");}
		else {out("Start at first SNP with full SNP window\n");}

		if(endSNP != 0) {out("End at SNP number: "); out(endSNP); out("\n");}
		else if(endSNPName != "") {out("End at SNP: "); out(endSNPName); out("\n");}
		else {out("End at last SNP with full SNP window\n");}

		if(jobNo != 0 && jobTotal != 0) {out("Job: "); out(jobNo); out(" of "); out(jobTotal); out("\n"); } 

		out("Window size: "); out(windowCMSize); out(" cM\n");
		out("Minimum number of SNPs in a window: "); out(windowMinSNPSize); out("\n");
		out("SNP step size: "); out(windowStepSize); out(" SNPs\n");
	};

	//if(missingCMZero) out("Treating centimorgan values set as 0 as missing\n");
	if(decreasingCM) out("Replacing decreasing centimorgan values with previous valid centimorgan value\n");

	//create analysis option and run analysis
	Analysis anAnalysis(filename, outputFilename, priorFilename, priorOutputFilename, plink, plinkPriorOptions, king, truffle, truffleOptions, gzip, ibdStitch, windowCMSize, missingCMZero, decreasingCM, windowMinSNPSize, windowStepSize, startSNP, endSNP, startSNPName, endSNPName, jobNo, jobTotal, noDomVar, priorOnly, keepTempFiles);

	if(!priorOnly) anAnalysis.runAnalysis();	

	time(&end);
	dif = difftime(end, start);
	out("\nRun time: "); out(getTime(dif)); out("\n\n");

	logFile.close();
};
