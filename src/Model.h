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


/*! \file Model.h
    \brief This file contains the models used for logistic regression.
    
*/

#ifndef __MODEL
#define __MODEL

#include "Data.h"

#include <map>

using namespace std;

#include "Data.h"

//! General Class for Models.
class Model
{
protected:

	map<unsigned int, double> parameters;
	
	SNPData * analysisSNP;		
	list<unsigned char> caseControls;

public:

	Model() : parameters() {};

	//setup parameter values and initValues of variables, to be over written in subclass of each model
	//where the variables created and added to the vector for the variables
	Model(map<unsigned int, double> paras) : parameters(paras) {};
	
	virtual ~Model() {};

	double getParameter(const unsigned int & no) const {return parameters.find(no)->second;};
	map<unsigned int, double> getParameters() const {return parameters;};
	void setNewParameters(map<unsigned int, double> & paras);
	void setAnalysisData(SNPData * anal) {analysisSNP = anal;};	
	void setCaseControls(const list<unsigned char> & ccs) {caseControls = ccs;};

	virtual void setParameter(const unsigned int & no, const double & value) {parameters[no] = value;};	
	virtual double negLogLikelihood() {return 0;};
	virtual void getGradientVector(map<unsigned int, double> & gradientVector, const bool & fitVarA, const bool & fitVarD) const {};
	virtual void getHessianMatrix(map<unsigned int, map<unsigned int, double> > & hessianMatrix, const bool & fitVarA, const bool & fitVarD) const {};
	virtual double getDeriv(unsigned int & paraToFit) const {return 0;};
	virtual double get2ndDeriv(unsigned int & paraToFit) const {return 1;};

};

//! Model for logistic regression on analysis SNP only.
class ModelIBD : public Model
{
private:

	double varA, varD; //additive and dominate variance parameters, actually varA/K^2 and varD/K^2
	IBDEstimates * priors; 
	IBDEstimates * posteriors;
	
public:

	ModelIBD() {};
	
	~ModelIBD() {};

	void setPriorsAndPosteriors(IBDEstimates * pri, IBDEstimates * post) {priors = pri; posteriors = post;};
	void setParameter(const unsigned int & no, const double & value);	
	double negLogLikelihood();
	
	void getGradientVector(map<unsigned int, double> & gradientVector, const bool & fitVarA, const bool & fitVarD) const;
	void getHessianMatrix(map<unsigned int, map<unsigned int, double> > & hessianMatrix, const bool & fitVarA, const bool & fitVarD) const;
	double getDeriv(unsigned int & paraToFit) const;
	double get2ndDeriv(unsigned int & paraToFit) const;
};

#endif
