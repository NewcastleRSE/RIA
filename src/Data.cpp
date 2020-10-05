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


/*! \file Data.cpp
    \brief This file contains the source for manipulating SNP data.
    
*/

#include <string>
#include <map>
#include <set>
#include <list>
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <cstdlib>
#include <algorithm>

using namespace std; // initiates the "std" or "standard" namespace

#include "main.h"
#include "Data.h"
#include "ModelFitting.h"
#include "Analysis.h"

//! Read in SNP data.
SNPData::SNPData(SNPWindow * snpWindow, const string & nm, const double & cm, const unsigned int & bp, const unsigned int & chr)
{
	unsigned int totalNoSubjects = snpWindow->getTotalNoSubjects();

	for(unsigned int subjectNo = 1; subjectNo <= totalNoSubjects; ++subjectNo)
	{
		noMinorAllelesAllSubjects.push_back(snpWindow->getNextNoOfMinorAlleles());
	};

	name = nm;
	cM = cm;
	bpPos = bp;	
	chromosome = chr;
};

//! Replaces SNP data of existing SNP object.
void SNPData::readInNewData(SNPWindow * snpWindow, const string & nm, const double & cm, const unsigned int & bp, const unsigned int & chr)
{	
	for(list<unsigned char>::iterator i = noMinorAllelesAllSubjects.begin(); i != noMinorAllelesAllSubjects.end(); ++i)
	{
		*i = snpWindow->getNextNoOfMinorAlleles();
	};

	name = nm;
	cM = cm;
	bpPos = bp;
	chromosome = chr;
};

//! Deletes all temp files.
void SNPWindow::deleteAllTemporaryFiles()
{
	analysis->deleteAllTemporaryFiles();
};

//! Gets the next number of minor alleles from file.
unsigned char SNPWindow::getNextNoOfMinorAlleles()
{
	int allele1, allele2;
	unsigned char noMinorAlleles = 0;

	//read in the next piece of data
	if(bitCount == 9)
	{
		
		readSNPData.read(buffer, 1);
		if(readSNPData.eof())
		{			
			outErr("Error: binary SNP file (.bed) is incomplete!\n");	
			deleteAllTemporaryFiles();
			exit(1);
		};
			
		aBit = buffer[0];
			
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

bool isChromosome1to22(unsigned int & chromosome)
{
	if(chromosome >= 1 && chromosome <= 22)	return true;

	return false;
};

//! Converts an u integer to a string
string toString(unsigned int & i)
{
	ostringstream aStringStream;
	aStringStream << i;

	return aStringStream.str();
};

//! Converts an double to a string
string toString(double & d)
{
	ostringstream aStringStream;
	aStringStream << d;

	return aStringStream.str();
};

//! Sets up SNP data from the .bim file and check cM data.
void SNPWindow::setUpSNPDesciptionData(string & filename, string & startSNPName, string & endSNPName)
{
	//try and find the binary map file, .bim, and read in data
	unsigned int length = filename.length();
	string mapFilename = filename.substr(0, length-4) + ".bim";

	ifstream readMapFile;
	readMapFile.open(mapFilename.c_str());
	if(!readMapFile.is_open())
	{
		outErr("Cannot read binary map file: "); outErr(mapFilename); outErr("!\n");
		deleteAllTemporaryFiles();
		exit(1);
	};

	string snpIdentifier;
	string prevSnpIdentifier = "";
	unsigned int chromosome, basePairPosition;
	double geneticDistance;
	string alleleName1, alleleName2;
	unsigned int snpNo = 1;

	//things for checking centimorgan data
	double firstCM = 0;
	double chromoCMMax = 0;
	unsigned int prevChr = 0;
	string SNPInfo;
	list<string> decreasingCMSNPs;
	list<string> decreasingCMInfoSNPs;
	list<unsigned int> chrsNoData;

	//loop thro' subjects and store the cases
	do{
		
		readMapFile >> chromosome >> snpIdentifier >> geneticDistance >> basePairPosition >> alleleName1 >> alleleName2;
		
		if(snpIdentifier != prevSnpIdentifier)
		{
				
				if(chromosome != prevChr)
				{
					if(prevChr != 0 && chromoCMMax == firstCM) chrsNoData.push_back(prevChr);
					firstCM = geneticDistance;
					chromoCMMax = geneticDistance;
				};
				
				if(geneticDistance < chromoCMMax)
				{
					if(decreasingCM)
					{
						geneticDistance = chromoCMMax;
					}
					else
					{
						SNPInfo = toString(chromosome) + "\t" + snpIdentifier + "\t" + toString(geneticDistance) + "\t" + toString(basePairPosition);
						decreasingCMInfoSNPs.push_back(SNPInfo);
						decreasingCMSNPs.push_back(snpIdentifier);
					};
				}
				else if(geneticDistance > chromoCMMax) chromoCMMax = geneticDistance;
				
				if(startSNPName == snpIdentifier) startSNP = snpNo;
				if(endSNPName == snpIdentifier) endSNP = snpNo;

				chromosomes.push_back(chromosome);				
				geneticDistances.push_back(geneticDistance);
				basePairs.push_back(basePairPosition);
				namesOfSNPs.push_back(snpIdentifier);

				prevChr = chromosome;
				++snpNo;
		};
		
		prevSnpIdentifier = snpIdentifier;
		
		if(!isChromosome1to22(chromosome))
		{
			outErr("RIA is designed for use with chromosomes 1 to 22 only!\n\n");
			outErr("Please remove all chromosome "); outErr(chromosome); outErr(" SNPs from your input files!\n\n");
			deleteAllTemporaryFiles();
			exit(1);
		};

	}while(!readMapFile.eof());

	readMapFile.close();

	if(prevChr != 0 && chromoCMMax == firstCM) chrsNoData.push_back(prevChr);

	out("Number of SNPs: "); out(geneticDistances.size()); out("\n");

	//output cM check info
	if(chrsNoData.size() > 0)
	{
		outErr("\nChromosome "); if(chrsNoData.size() > 1) out("s ");

		list<unsigned int>::const_iterator cnd = chrsNoData.begin();
		unsigned int chrNo = 1;
		while(cnd != chrsNoData.end())
		{
			outErr(*cnd);
			
			if(chrNo == (chrsNoData.size() - 1)) outErr(" and ");
			else if(chrNo < chrsNoData.size()) outErr(", ");

			++cnd;
			++chrNo;			
		};
		if(chrsNoData.size() > 1) outErr(" do not have any centimorgan data!\n");
		else outErr(" does not have any centimorgan data!\n");

		deleteAllTemporaryFiles();
		exit(1);
	}
	else if(decreasingCMSNPs.size() > 0)
	{
		outErr("\n");
		
		if(decreasingCMSNPs.size() == 1)
		{
			outErr("There is one SNP with a decreasing centimorgan value!\n");
		}
		else
		{
			outErr("There are "); outErr(decreasingCMSNPs.size()); outErr(" SNPs with a decreasing centimorgan value!\n");
		};
		
		outErr("\n");
	
		outErr("SNPs with a decreasing centimorgan value are listed below:\n");
		unsigned int count = 1;
		for(list<string>::const_iterator dcs = decreasingCMInfoSNPs.begin(); dcs != decreasingCMInfoSNPs.end(); ++dcs)
		{
			outErr(*dcs); outErr("\n");
			if(count > 50) {outErr("and so on...\n"); break;};
			count++;
		};			
		
		string filenameD = "ria-decreasing-cm-snps.txt";
		outErr("\n");
		outErr("These SNPs are listed in the file "); outErr(filenameD); outErr("\n");
			
		ofstream probSNPs(filenameD.c_str());

		for(list<string>::const_iterator dcs = decreasingCMSNPs.begin(); dcs != decreasingCMSNPs.end(); ++dcs)
		{
			probSNPs << *dcs << "\n";
		};
		probSNPs.close();

		outErr("\n");
		outErr("Please do one of the following:\n");
		outErr("1) Correct your centimorgan SNP data to be increasing for each chromosome.\n");
		outErr("2) Remove the problematic SNPs, for example:\n\t plink --noweb --bfile "); outErr(filename.substr(0, length-4)); outErr(" --exclude ria-decreasing-cm-snps.txt");
			outErr(" --make-bed --out "); outErr(filename.substr(0, length-4)); outErr("-fixed\n");
		outErr("3) Use the \"-decreasing-cm\" option to use the previous valid centimorgan value for SNPs with a decreasing centimorgan value.\n\n");	

		deleteAllTemporaryFiles();

		exit(1);
	};
};

//! Sets up family case/control data from the .fam file.
void SNPWindow::setupCaseControls(string & fname, const bool & makePostFamFile)
{
	//try and find the family file and read in data
	unsigned int length = fname.length();
	string famFilename = fname.substr(0, length-4) + ".fam";

	ifstream readFamilyFile;
	readFamilyFile.open(famFilename.c_str());
	if(!readFamilyFile.is_open())
	{
		outErr("Cannot read family file: "); outErr(famFilename); outErr("!\n");
		deleteAllTemporaryFiles();
		exit(1);
	};

	ofstream postOutFamFile;
	if(makePostFamFile)
	{
		postOutFamFile.open(postFamFile.c_str());
	};

	string famID, indivID, FatherId, MotherID, sexID, famIndivID;
	string prevFamIndivID = "";
	int phenoType;
	unsigned int noCases = 0;

	//loop thro' subjects and store the cases
	do{		
		readFamilyFile >> famID >> indivID >> FatherId >> MotherID >> sexID >> phenoType;
		famIndivID = famID + "-" + indivID;

		//do not duplicate the last row
		if(famIndivID != prevFamIndivID) 
		{
			if(phenoType == 2)
			{
					noCases++;
					caseControls.push_back(2);
					if(makePostFamFile) postOutFamFile << famID << " " << indivID << " " << FatherId << " " << MotherID << " " << sexID << " " << phenoType << "\n";
			}			
			else if(phenoType == 1) caseControls.push_back(1);
			else caseControls.push_back(0);
		};

		prevFamIndivID = famIndivID;
	}while(!readFamilyFile.eof());

	if(makePostFamFile) postOutFamFile.close();
	readFamilyFile.close();

	totalNoSubjects = caseControls.size(); //needed for reading in data later

	out("Number of cases: "); out(noCases); out("\n");
	if(noCases != totalNoSubjects) {out("Number of unused controls: "); out((totalNoSubjects-noCases)); out("\n");};

};

//! Writes binary data.
void writeBedGenotype(ofstream & psBedFile, unsigned int & psBitCount, int & aBit, const bool & allele1, const bool & allele2)
{

	aBit = aBit >> 1; //shift bits to the right
	if(allele1) aBit = aBit | 128; 
	aBit = aBit >> 1;  //shift bits to the right
	if(allele2) aBit = aBit | 128;
	
	psBitCount += 2;

	//write to file if byte is finished
	if(psBitCount == 8)
	{
		//write to file
		char buffer[1];
		buffer[0] = aBit;
		psBedFile.write(buffer, 1);
		psBitCount = 0;
		aBit = 0;
	};

};

//! Writes last byte and resets the Bit and bit counter
void writeLastByteBeforeNextSNP(ofstream & psBedFile, unsigned int & psBitCount, int & aBit)
{
	if(psBitCount == 0) return; //last byte may already be written

	//shift right bits over
	while(psBitCount < 8)
	{
		aBit = aBit >> 1; //shift bits to the right
		psBitCount++;
	};
	
	//write to file
	char buffer[1];
	buffer[0] = aBit;
	psBedFile.write(buffer, 1);
	psBitCount = 0;
	aBit = 0;
};

//! Creates files used by KING to calculate posterior IBD estimates.
void SNPWindow::createFilesForPosteriorCalc(string & outputFilename)
{
	createFamFileForPosteriorCalc(outputFilename);
	createBimFileForPosteriorCalc(outputFilename);

	unsigned int psBitCount = 0;
	int aBit = 0;
	bool allele1;
	bool allele2;

	string outputBedFilename = outputFilename + ".bed";
	ofstream outBedFile(outputBedFilename.c_str(), ios::binary);

	//write out initial binary pedigree file bytes, first 2 byte are magic numbers the third to indicate SNP major (subjects x SNPs)
	char buffer[3];	
	buffer[0] = 108;
	buffer[1] = 27;
	buffer[2] = 1;
	outBedFile.write(buffer, 3);

	list<unsigned char>::const_iterator cc; //case/control status, 1=case
	
	//loop through SNPs
	for(map<unsigned int, SNPData *>::const_iterator snp = window.begin(); snp != window.end(); ++snp)
	{	
		cc = caseControls.begin();

		//loop thro subjects for this SNP
		for(list<unsigned char>::const_iterator i = snp->second->noMinorAllelesAllSubjects.begin(); i != snp->second->noMinorAllelesAllSubjects.end() && cc != caseControls.end(); ++i)
		{
			//if subject is a case record the SNP info
			if(*cc == 2)
			{
				if(*i == 0)
				{
					allele1 = true;
					allele2 = true;
				}
				else if(*i == 1)
				{
					allele1 = false;
					allele2 = true;
				}
				else if(*i == 2)
				{
					allele1 = false;
					allele2 = false;
				}
				else if(*i == 3) //missing data
				{
					allele1 = true;
					allele2 = false;
				};

				writeBedGenotype(outBedFile, psBitCount, aBit, allele1, allele2);
			};

			++cc;
		};

		//start new byte as starting a new SNP
		writeLastByteBeforeNextSNP(outBedFile, psBitCount, aBit);
	};

	outBedFile.close();
};

//! Creates .fam file used by KING to calculate posterior IBD estimates.
void SNPWindow::createFamFileForPosteriorCalc(string & outputFilename)
{
	//make copy of .fam file for use with posterior calculations	 
	string copyFamCommand = cpStr + postFamFile + " " + outputFilename + ".fam " + endCommand;
	unsigned int cpErrorCode = system(copyFamCommand.c_str());
	if(cpErrorCode != 0)
	{
		outErr("Problem copying (error code "); outErr(cpErrorCode); outErr(") .fam file for posterior calculations with:\n\n ");
		outErr(copyFamCommand); outErr("\n\n");
		exit(cpErrorCode);		
	};

};

//! Creates .bim file used by KING to calculate posterior IBD estimates.
void SNPWindow::createBimFileForPosteriorCalc(string & outputFilename)
{
	string bimOutputFilename = outputFilename + ".bim";
	ofstream bimFile(bimOutputFilename.c_str());

	//unsigned int basePairPosition = 10;

	//Loop thro' SNPs
	for(map<unsigned int, SNPData *>::const_iterator snp = window.begin(); snp != window.end(); ++snp)
	{
		bimFile << snp->second->chromosome << " " << snp->second->name << " " << snp->second->cM << " " << snp->second->bpPos << " A G\n";
			
		//basePairPosition += 10;
	};

	bimFile.close();
};

//! Adds the data for the next SNP to the SNP window.
void SNPWindow::addNextSNPDataToWindow()
{
	//create object to store data for the one SNP
	bitCount = 9; //move reading data onto the next byte
	SNPData * someSNPData;
	
	someSNPData = getNewSNPData();
	
	if(!newSNPObject)
	{
		someSNPData->readInNewData(this, *nameSNPiter, *cmSNPiter, *bpSNPiter, *chrSNPiter);
		nameSNPiter++;
		cmSNPiter++;
		bpSNPiter++;
		chrSNPiter++;
	};

	window[addingSNPNo++] = someSNPData;
};

//! Creates the SNP window.
SNPWindow::SNPWindow(string & filename, double & wcMs, unsigned int & wmssnps, unsigned int & wstps, bool & misscm, bool decm, unsigned int & stSNPno, unsigned int & eSNPno, string & stSNPname, string & eSNPname, string & cp, string & rm, string & end, string & pff, const bool & priorsGiven, bool & priorOnly, unsigned int & jobNo, unsigned int & jobTotal, Analysis * an)
	: geneticDistances(), window(), windowSizes(), windowCMSize(wcMs), missingCMZero(misscm), decreasingCM(decm), windowStepSize(wstps), startSNP(stSNPno), endSNP(eSNPno), analysisSNPNo(1), addingSNPNo(1), newChr(false), validWindow(true), firstCMChr(0), cpStr(cp), rmStr(rm), endCommand(end), postFamFile(pff), analysis(an), startSNPJobNo(0), endSNPJobNo(0), noCalcs(1)
{	

	setupCaseControls(filename, priorsGiven);

	//used for reading in SNP binary data
	one = '\1';
	bitCount = 9;
	
	setUpSNPDesciptionData(filename, stSNPname, eSNPname); 

	if(!priorOnly) 
	{
		setWindowMinSNPSize(wmssnps);

		if(jobNo != 0) setStartAndEndSNPBasedOnJobNo(jobNo, jobTotal);

		setUpFirstWindow(filename);
	};

	out("\n");
};

//! Moves through SNP data by one SNP.
void SNPWindow::advanceSNPData()
{
	//a new byte is started after each SNP, 4 subject genotypes per byte,
	// so the no. of bytes is rounded up when divided by 4
	unsigned int bufferSize;
	if(totalNoSubjects%4 == 0) bufferSize = totalNoSubjects/4;
	else bufferSize = (unsigned int)(totalNoSubjects/4) + 1;

	char buffer[1];
	for(unsigned int i = 1; i <= bufferSize; ++i)
	{
		readSNPData.read(buffer, 1);
	};

};

//! Estimate total posterior calculations to be made.
void SNPWindow::estimateTotalPostCalcs(const unsigned int & jobTotal)
{
	list<double>::const_iterator gd = geneticDistances.begin(); 
	list<unsigned int>::const_iterator chr = chromosomes.begin();

	double windowCMSizeHalf = windowCMSize*0.5; 
	double startGD = 0;
	unsigned int prevChr = 999;
	double prevGD = 0;
	unsigned int snpNo = 1;
	unsigned int removeSNPs = 0;
	bool firstChr = true;
	double lastOKGD = 0;

	list<double> lastOKGDs;
	
	//get genetic distance of last valid SNP in each chr
	while(gd != geneticDistances.end() && chr != chromosomes.end())
	{
		prevChr = *chr;
		prevGD = *gd;
		++gd; ++chr; 

		if(gd == geneticDistances.end() || (chr != chromosomes.end() && *chr != prevChr))
		{	
			lastOKGD = prevGD - windowCMSizeHalf;
			lastOKGDs.push_back(lastOKGD);			
		};

	};

	prevChr = 999;
	prevGD = 0;
	gd = geneticDistances.begin(); 
	chr = chromosomes.begin();
	list<double>::const_iterator logd = lastOKGDs.begin();

	if(logd != lastOKGDs.end()) lastOKGD = *logd;

	while(gd != geneticDistances.end() && chr != chromosomes.end() 
		&& ((endSNP == 0 && endSNPJobNo == 0) || (endSNP != 0 && snpNo <= endSNP) || (endSNPJobNo != 0 && snpNo <= endSNPJobNo)))
	{
		if(*chr != prevChr)
		{
			startGD = *gd;
			if(!firstChr)
			{
				++logd;
				if(logd != lastOKGDs.end()) lastOKGD = *logd;
			};

			firstChr = false;
		};

		if(!(
			*gd >= (startGD + windowCMSizeHalf)
		    && *gd <= lastOKGD
			&& ((startSNP == 0 && startSNPJobNo == 0)
				|| (startSNPJobNo == 0 && startSNP != 0 && snpNo >= startSNP) 
				|| (startSNPJobNo != 0 && snpNo >= startSNPJobNo))			
			))
		{
			removeSNPs++; //remove at start chr				
		};

		prevChr = *chr;
		prevGD = *gd;
		++gd; ++chr; ++snpNo;
	};

	noCalcs = (snpNo - 1 - removeSNPs)/(double)(windowStepSize);

	if(noCalcs < 1) noCalcs = 1;
};

//! Outputs percentage progress to screen.
void SNPWindow::outputProgress(const unsigned int & calcCount)
{
	
	double percent = (double)calcCount/noCalcs;

	percent = ((double)((int)(percent*10000 + 0.5))/100.0);

	if(percent > 100) percent = 100;
	
	screenOut("\r                             "); screenFlush("");
	screenOut("\rProgress: "); screenOut(percent); screenOut("%"); screenFlush("");

};

//! Test code for outputting window info.
void SNPWindow::displayWindowInfo()
{
	
	cout << analysisSNPNo << " " << analysisSNP->cM <<" " << analysisSNP->bpPos <<" " << analysisSNP->name << ": ";
	map<unsigned int, SNPData *>::const_iterator start = window.begin();
	map<unsigned int, SNPData *>::const_reverse_iterator end = window.rbegin();

	end++; //end SNP is not in the window

	double startCM = 0;
	double endCM = 0;

	if(start != window.end())
	{
		cout << start->first << " " << start->second->name << " ";
		startCM = start->second->cM;
		cout << startCM << " ";
	};

	if(end != window.rend())
	{
		cout << "| " << end->first << " " << end->second->name << " ";
		endCM = end->second->cM;
		cout << endCM << " ";
	};

	if(start != window.end() && end != window.rend()) cout << " | size = " << endCM - startCM << ", SNPs = " << window.size();
	cout << "\n";
};

//! Gets first SNP in initial window.
unsigned int SNPWindow::getStartCMWindowSNP(unsigned int & startAnalysisSNPno)
{
	//get cM position of start SNP
	double posStartCM;
	unsigned int snpNo = 1;
	unsigned int startChr = 0;

	double halfWindowCMSizePlus = windowCMSize*0.5 + *geneticDistances.begin();

	bool updateStartAnalSNP = (startAnalysisSNPno == 0);
	if(updateStartAnalSNP) updateStartAnalSNP = 1;

	list<unsigned int>::const_iterator chr = chromosomes.begin();
	for(list<double>::const_iterator i = geneticDistances.begin(); i != geneticDistances.end() && chr != chromosomes.end(); ++i, ++chr)
	{
		if(!updateStartAnalSNP && startAnalysisSNPno == snpNo)
		{
			posStartCM = *i;
			startChr = *chr;
			break;
		}
		else if(updateStartAnalSNP && *i >= halfWindowCMSizePlus)
		{
			posStartCM = *i;
			startChr = *chr;
			startAnalysisSNPno = snpNo;			
			break;
		};

		snpNo++;
	};
	
	double startCMWindow; //LHS SNP window cM start position
	
	if(posStartCM < halfWindowCMSizePlus)
	{
		outErr("A full window of size "); outErr(windowCMSize); outErr(" cM is not possible with the given start SNP!\n");
		outErr("Please set the start SNP to a SNP with at least "); outErr(windowCMSize*0.5); outErr(" cM from the first SNP in a chromosome, or set a smaller SNP window size\n");			
		outErr("Alternatively, do not set the start SNP to automatically choose the first possible window.\n");
		outErr("Also check that the centimorgan SNP data is complete and correct!\n");
		deleteAllTemporaryFiles();
		exit(1);
	}
	else startCMWindow = posStartCM - windowCMSize*0.5;

	chr = chromosomes.begin();
	snpNo = 1;
	for(list<double>::const_iterator i = geneticDistances.begin(); i != geneticDistances.end() && chr != chromosomes.end(); ++i, ++chr)
	{
		if(*i >= startCMWindow && *chr == startChr) return snpNo;
		snpNo++;
	};

	outErr("Problem setting up first cM analysis window!\n");
	outErr("Please check that the centimorgan SNP data is set correctly.\n");
	outErr("Also that the start SNP is set correctly if set. Do not set to automatically choose the first possible window.\n");
	outErr("Also check that the centimorgan SNP data is complete and correct!\n");

	deleteAllTemporaryFiles();
	exit(1);

	return 1;//this should not be possible
};

//! Moves SNP window onto the chosen first SNP, but does not read data in for first SNP in window.
void SNPWindow::advanceToFirstWindow(unsigned int & analysisSNPNo)
{
	unsigned int snpCount = 1;

	unsigned int snpToStopAt = getStartCMWindowSNP(analysisSNPNo);

	do{
		if(snpCount == snpToStopAt) break;

		//read SNP data in for the previous SNP that will not be used
		advanceSNPData();
		addingSNPNo++; //keep this in sync as if the SNP data were added

		snpCount++;

		nameSNPiter++;
		cmSNPiter++;
		bpSNPiter++;
		chrSNPiter++;
		analysisCMSNP++;
		analysisChrSNP++;
		endCMSNP++; endChr++;
	}while(cmSNPiter !=  geneticDistances.end());


	unsigned int noTimesToAdvAnalysisCM = analysisSNPNo - snpToStopAt;

	for(unsigned int i = 1; i <= noTimesToAdvAnalysisCM; ++i) {analysisCMSNP++; analysisChrSNP++;};

};

//! Tests end of SNP window are in the same chr.
bool SNPWindow::isWindowInOneChr()
{
	unsigned int chrStart = 0;
	unsigned int chrEnd = 0;
	if(window.begin() != window.end())
	{
			chrStart = window.begin()->second->chromosome;
			map<unsigned int, SNPData *>::const_reverse_iterator ec = window.rbegin();
			if(ec != window.rend()) chrEnd = ec->second->chromosome;
	};

	return chrStart == chrEnd;
};

//! Tests end of SNP window are in the same chr.
unsigned int SNPWindow::getLastWindowChr()
{
	unsigned int chrEnd = 0;
	
	map<unsigned int, SNPData *>::const_reverse_iterator ec = window.rbegin();
	if(ec != window.rend()) chrEnd = ec->second->chromosome;
	
	return chrEnd; 
}

//! Update the LHS of the SNP window based on BP position.
void SNPWindow::updateStartSNP()
{		
	double startWindowCM = 0;
	if(*analysisCMSNP > windowCMSize*0.5 && !newChr) startWindowCM = *analysisCMSNP - windowCMSize*0.5;	
	
	unsigned int nextChr = 0;
	if(newChr)
	{
		nextChr = getLastWindowChr();
	};

	//remove all SNPs in window	that are outside the new LHS cM SNP window boundary
	map<unsigned int, SNPData *>::iterator aSNP = window.begin();

	while(aSNP != window.end() && ((!newChr && aSNP->second->cM < startWindowCM) || (newChr && aSNP->second->chromosome != nextChr)))
	{
		spareSNPs.push(aSNP->second);
		window.erase(aSNP);
		aSNP = window.begin();
	};	
};

//! Add SNPs to the required window based on CM size until all the RHS is filled with SNP data. 
bool SNPWindow::updateEndSNP()
{	
	bool isOK = true;
	validWindow = true;

	double endWindowCM = *analysisCMSNP + windowCMSize*0.5;	
	double prev = *endCMSNP;
	double prevEndChr = *endChr;

	//run iterator thro' SNPs until we go past the end of the window or come to the end of the file
	do{
		
		if(endCMSNP == geneticDistances.end() || endChr == chromosomes.end()) return false;
		else if(*endCMSNP > endWindowCM && *endChr == *analysisChrSNP) break; //window big enough and not switched chrs		
		else if(endCMSNP != geneticDistances.end() && endChr != chromosomes.end())
		{
			prev = *endCMSNP;
			prevEndChr = *endChr;
			endCMSNP++;
			endChr++;		
		};

		if(endCMSNP != geneticDistances.end() && endChr != chromosomes.end())
		{
			addNextSNPDataToWindow();

			if(prevEndChr != *endChr) //new chr
			{
				validWindow = false; 
				firstCMChr = *endCMSNP;				
				break;
			};
		};
	}while(true);

	return isOK;
};

//! Moves SNP window onto first SNP in the next chr.
void SNPWindow::advanceToFirstWindowInChr()
{
	unsigned int chrEnd = 0;
	if(window.begin() != window.end())
	{
		map<unsigned int, SNPData *>::const_reverse_iterator ec = window.rbegin();
		if(ec != window.rend()) chrEnd = ec->second->chromosome;
	};

	if(chrEnd == 0)
	{
		validWindow = false;
		return;
	};

	list<unsigned int>::const_iterator c = chromosomes.begin();
	unsigned int firstSNPNoInNextChr = 1;

	while(c != chromosomes.end() && *c != chrEnd)
	{
		firstSNPNoInNextChr++;
		c++;
	};

	if(c == chromosomes.end() || firstSNPNoInNextChr < analysisSNPNo)
	{
		validWindow = false;
		return;
	};

	double halfWindowSizePlus = windowCMSize*0.5 + firstCMChr;


	while(analysisCMSNP != geneticDistances.end() && (analysisSNPNo < firstSNPNoInNextChr || *analysisCMSNP < halfWindowSizePlus))
	{
		analysisSNPNo++;	

		//move the window along one SNP in the list of all SNPs
		analysisCMSNP++;
		analysisChrSNP++;
	};

	analysisSNPNo--; //move back one
	analysisCMSNP--;
	analysisChrSNP--;

	if(analysisCMSNP == geneticDistances.end()) validWindow = false;
};

//! Sets up the SNP IDs for the start and end SNPs of the first SNP window.
void SNPWindow::setUpStartAndEndSNPs()
{
	if(startSNP > geneticDistances.size())
	{
		outErr("The start SNP must not be greater than the total number of SNPs!\n");
		deleteAllTemporaryFiles();
		exit(1);
	};

	analysisCMSNP = geneticDistances.begin();
	analysisChrSNP = chromosomes.begin();

	//Set SNP to start analysis
	analysisSNPNo = startSNP; 

	endCMSNP = geneticDistances.begin();
	endChr = chromosomes.begin();
	nameSNPiter = namesOfSNPs.begin();
	cmSNPiter = geneticDistances.begin();
	bpSNPiter = basePairs.begin();
	chrSNPiter = chromosomes.begin();

	advanceToFirstWindow(analysisSNPNo);

	//add SNP data for the first SNP in window
	addNextSNPDataToWindow();
	
	//move to end SNP according to the cM window size and add data to the window
	if(!updateEndSNP() || !validWindow)
	{
		outErr("Unable to set up even one analysis SNP window!\n");
		outErr("Check the start SNP is not too close to the end of a chromosome, and that the cM window size is not too big.\n");
		outErr("Also check that the centimorgan SNP data is complete and correct!\n");
		deleteAllTemporaryFiles();
		exit(1);
	};

	map<unsigned int, SNPData *>::const_iterator ws = window.find(analysisSNPNo);
	if(ws != window.end()) analysisSNP = ws->second;
	else
	{
		outErr("Problem setting up initial SNP window, please check cM data!\n");
		deleteAllTemporaryFiles();
		exit(1);
	};
};

//! Sets up the start and end SNP nos based on job n of N type things.
void SNPWindow::setStartAndEndSNPBasedOnJobNo(unsigned int & jobNo, unsigned int & jobTotal)
{
	double noSNPs = (double)(basePairs.size());

	double stepSize = noSNPs/(double)(jobTotal);

	//check job options are acceptable
	if(stepSize <= (double)(windowStepSize))
	{
		outErr("The SNP window step size, "); outErr(windowStepSize); outErr(", is not allowed to be greater than the number of SNPs spanned in one job, ");
		outErr("which is about "); outErr(stepSize); outErr(" SNPs!\n");

		deleteAllTemporaryFiles();
		exit(1);
	};

	//round up start SNP to be multiple of windowStepSize, then plus 1
	startSNPJobNo = (unsigned int)(stepSize*(jobNo-1) + 1);

	endSNPJobNo = (unsigned int)(stepSize*jobNo);

};

bool SNPWindow::isAnalysisSNPInJobWindow()
{
	return startSNPJobNo == 0 || (analysisSNPNo >= startSNPJobNo && analysisSNPNo <= endSNPJobNo);
};

//! Updates start and end SNP in the analysis object.
void SNPWindow::updateStartAndEndSNP(unsigned int & startSNPAnal, unsigned int & endSNPAnal)
{
	startSNPAnal = startSNP;
	endSNPAnal = endSNP;
};

//! Creates the first SNP window, opens SNP file and reads in data for the first window.
void SNPWindow::setUpFirstWindow(string & filename)
{
	//try and find the binary pedigree file, .bed, and read in data for the first window
	readSNPData.open(filename.c_str(), ios::binary);
	
	if(!readSNPData.is_open())
	{
		outErr("Cannot read binary pedigree file: "); outErr(filename); outErr("!\n");
		deleteAllTemporaryFiles();
		exit(1);
	};

	char buffer[3];
	readSNPData.read(buffer, 3);

	//check the plink magic numbers for the file type
	//3rd number indicates format of genotype data, 1 => subjects x SNPs, 0 => SNPs x subjects
	unsigned int magicNumber1 = buffer[0];
	unsigned int magicNumber2 = buffer[1];

	if(magicNumber1 != 108 || magicNumber2 != 27)
	{
		outErr( "Detected an old version .bed file!\n");
		outErr("Please use PLINK to update the .bed file.\n");
			
		readSNPData.close();
		deleteAllTemporaryFiles();
		exit(1);
	};

	//determine binary file type
	unsigned int mode = buffer[2];
	if(mode == 0)
	{
		outErr("The binary pedigree file must be in SNP-major mode!\n");
		outErr("Please use PLINK to update the .bed file.\n");
			
		readSNPData.close();	
		deleteAllTemporaryFiles();
		exit(1);
	};

	setUpStartAndEndSNPs();

};

//! Gets new SNP data from spare SNP data objects or creates a new SNP Data object.
SNPData * SNPWindow::getNewSNPData()
{
	if(spareSNPs.size() == 0)
	{	
		newSNPObject = true;
		return new SNPData(this, *(nameSNPiter++), *(cmSNPiter++), *(bpSNPiter++), *(chrSNPiter++));
	};

	SNPData * snpData = spareSNPs.top();
	spareSNPs.pop();
	newSNPObject = false;

	return snpData;
};

//! Moves the window on for the next SNP.
bool SNPWindow::nextWindow()
{
	bool atEndOfData = false;

	validWindow = true;

	//setup new analysis SNP
	map<unsigned int, SNPData *>::const_iterator ws = window.find(analysisSNPNo);
	if(ws != window.end()) analysisSNP = ws->second;
	else
	{
		outErr("Problem setting up SNP window, please check cM data!\n");
		deleteAllTemporaryFiles();
		exit(1);
	};
	

	if(newChr)
	{
		advanceToFirstWindowInChr();
	}
	else
	{
		unsigned int prevAnalChr = *analysisChrSNP;
		//move the window along by the number of SNPs in given number step size in the list of all SNPs
		for(unsigned int i = 1; i <= windowStepSize; ++i)
		{
			//setup new analysis SNP number
			analysisSNPNo++;

			//move the window along one SNP in the list of all SNPs
			analysisCMSNP++;
			analysisChrSNP++;			
			if(analysisCMSNP == geneticDistances.end()) return true; //at end of data nothing more to do
			if(prevAnalChr != *analysisChrSNP) break;
		};

		if(prevAnalChr != *analysisChrSNP)
		{			
			firstCMChr = *analysisCMSNP;	
			advanceToFirstWindowInChr();
		};
	};

	
	if(!updateEndSNP()) return true; //return true if at end of the data
	
	//setup new analysis SNP	
	ws = window.find(analysisSNPNo);
	if(ws != window.end()) analysisSNP = ws->second;
	else
	{
		outErr("Problem creating the SNP window, please check cM data!\n");
		deleteAllTemporaryFiles();
		exit(1);
	};

	updateStartSNP(); //set LHS snip window based on cM value of analysis SNP
	
	if(!isWindowInOneChr()) //at end of chr
	{
		validWindow = false;
		newChr = true;
	}
	else
	{
		newChr = false;
	};

	if(endSNP != 0 && analysisSNPNo > endSNP) atEndOfData = true; //no more analysis to do

	return atEndOfData;
};

bool myfunction (unsigned int i, unsigned int j) { return (i<j); }

//! Displays stats of the window size, mean, sd etc.
void SNPWindow::displayWindowStats()
{
	double winMean = 0;
	double winMedian = 0;
	unsigned int halfNo = windowSizes.size()/2;
	if(halfNo == 0) halfNo = 1;
	double winSD = 0;
	double minWin = (*(windowSizes.begin()));
	double maxWin = (*(windowSizes.begin()));

	multiset<unsigned int> orderedSizes;

	//sort(windowSizes.begin(), windowSizes.end(), myfunction);

	bool even = windowSizes.size()%2 == 0;
	unsigned int count = 1;

	for(list<unsigned int>::const_iterator ws = windowSizes.begin(); ws != windowSizes.end(); ++ws)
	{	
		orderedSizes.insert(*ws);
		winMean += *ws;
		if(*ws < minWin) minWin = *ws;
		if(*ws > maxWin) maxWin = *ws;
	};

	winMean /= (double)(windowSizes.size());
	
	for(multiset<unsigned int>::const_iterator os = orderedSizes.begin(); os != orderedSizes.end(); ++os)
	{
		if(count == halfNo)
		{
			winMedian = *os;
		}
		else if(count == (halfNo + 1) && even) winMedian = (winMedian + *os)*0.5;
		else if(count > (halfNo + 1)) break;

		count++;
	};


	for(list<unsigned int>::const_iterator ws = windowSizes.begin(); ws != windowSizes.end(); ++ws)
	{		
		winSD += (*ws - winMean)*(*ws - winMean);		
	};

	winSD /= (double)(windowSizes.size());
	winSD = sqrt(winSD);

	out("\n");
	out(windowCMSize); out(" cM windows - number of SNPs summary statistics:\n");
	if(windowSizes.size() != 0)
	{
		out("Mean: "); out(winMean); out("\n");
		out("Median: "); out(winMedian); out("\n");
		out("Standard deviation: "); out(winSD); out("\n");
		out("Range: ("); out(minWin); out(", "); out(maxWin); out(")\n\n");
	} else {
		out("No results were calculated for this analysis.");
	};

};
