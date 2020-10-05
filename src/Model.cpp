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


/*! \file Model.cpp
    \brief This file contains the source of models used for logistic regression.
    
*/

#include <map>
#include <iostream>
#include <math.h>

using namespace std; // initiates the "std" or "standard" namespace

#include "Model.h"
#include "Data.h"
#include "ModelFitting.h"

//! Sets the values of all of the parameters.
void Model::setNewParameters(map<unsigned int, double> & paras)
{
	for(map<unsigned int, double>::const_iterator p = paras.begin(); p != paras.end(); ++p)
	{
		setParameter(p->first, p->second);
	};
};

//! Sets parameters.
void ModelIBD::setParameter(const unsigned int & no, const double & value)
{
	parameters[no] = value;

	switch(no)
	{
		case 1:
			varA = value;
			return;
		case 2:
			varD = value;
			return;
	};
};

//! Returns negative log (base 10!) likelihood for single-SNP log regression.
double ModelIBD::negLogLikelihood()
{
	double ans = 0;
	double priorBit, postBit;
	bool isSomeData = false;

	//loop through prior and posterior estimates and calculate neglog like
	map<string, IBDData *>::const_iterator pri = priors->snpIBDEstimates.begin();
	map<string, IBDData *>::const_iterator post = posteriors->snpIBDEstimates.begin();

	do{
		//check prior and posteriors are for the same pair
		if(post == posteriors->snpIBDEstimates.end())
		{
			//find correct pair in posteriors
			post = posteriors->snpIBDEstimates.find(pri->first);
		}
		else if(post != posteriors->snpIBDEstimates.end() && pri->first != post->first)
		{
			//find correct pair in posteriors
			post = posteriors->snpIBDEstimates.find(pri->first);			
		};

		if(post != posteriors->snpIBDEstimates.end())
		{
			//ignore pairs with bad estimates
			if(!pri->second->error && !post->second->error)
			{
				isSomeData = true;
				priorBit = 1 + pri->second->f2j*(varA + varD) + 0.5*pri->second->f1j*varA;
				postBit = 1 + post->second->f2j*(varA + varD) + 0.5*post->second->f1j*varA;
				ans += log10(priorBit) - log10(postBit);			
			};
		};

		++pri;
		if(post != posteriors->snpIBDEstimates.end()) ++post;
	}while(pri != priors->snpIBDEstimates.end());

	if(!isSomeData) return log10(ans); //returns a nan (=log10(0)) if there is no data so does not converge to such data, will give NA in results

	return ans;
};

//! Get gradient vector for 2 parameter model.
void ModelIBD::getGradientVector(map<unsigned int, double> & gradientVector, const bool & fitVarA, const bool & fitVarD)  const
{
	//fit varA and varD
	double ans[2] = {0, 0};

	double priorBit, postBit;
	double priNum, postNum;

	//loop through prior and posterior estimates and calculate neglog like
	map<string, IBDData *>::const_iterator pri = priors->snpIBDEstimates.begin();
	map<string, IBDData *>::const_iterator post = posteriors->snpIBDEstimates.begin();

	do{
		//check prior and posteriors are for the same pair
		if(post == posteriors->snpIBDEstimates.end())
		{
			//find correct pair in posteriors
			post = posteriors->snpIBDEstimates.find(pri->first);
		}
		else if(post != posteriors->snpIBDEstimates.end() && pri->first != post->first)
		{
			//find correct pair in posteriors
			post = posteriors->snpIBDEstimates.find(pri->first);			
		};

		if(post != posteriors->snpIBDEstimates.end())
		{
			//ignore pairs with bad estimates
			if(!pri->second->error && !post->second->error)
			{
				//dl/dV_A
				priNum = pri->second->f2j + 0.5*pri->second->f1j;
				postNum = post->second->f2j + 0.5*post->second->f1j;
				priorBit = priNum/(1 + priNum*varA + pri->second->f2j*varD);
				postBit = postNum/(1 + postNum*varA + post->second->f2j*varD);

				ans[0] += priorBit - postBit;

				//dl/dV_D				
				priorBit = pri->second->f2j/(1 + priNum*varA + pri->second->f2j*varD);
				postBit = post->second->f2j/(1 + postNum*varA + post->second->f2j*varD);

				ans[1] += priorBit - postBit;
			};
		};

		++pri;
		if(post != posteriors->snpIBDEstimates.end()) ++post;
	}while(pri != priors->snpIBDEstimates.end());
	
	gradientVector[1] = ans[0];	
	gradientVector[2] = ans[1];	

};

//! Returns the 2nd derivative w.r.t. chosen parameters of the negative log likelihood.
void ModelIBD::getHessianMatrix(map<unsigned int, map<unsigned int, double> > & hessianMatrix, const bool & fitVarA, const bool & fitVarD) const
{
	//Always fitting beta0 and beta1 for 1 SNP model
	double ans[2][2] = {{0, 0}, {0, 0}}; //col, row
	
	double priorBit, postBit;
	double priNum, postNum, priDenom, postDenom;

	//loop through prior and posterior estimates and calculate neglog like
	map<string, IBDData *>::const_iterator pri = priors->snpIBDEstimates.begin();
	map<string, IBDData *>::const_iterator post = posteriors->snpIBDEstimates.begin();

	do{
		//check prior and posteriors are for the same pair
		if(post == posteriors->snpIBDEstimates.end())
		{
			//find correct pair in posteriors
			post = posteriors->snpIBDEstimates.find(pri->first);
		}
		else if(post != posteriors->snpIBDEstimates.end() && pri->first != post->first)
		{
			//find correct pair in posteriors
			post = posteriors->snpIBDEstimates.find(pri->first);			
		};

		if(post != posteriors->snpIBDEstimates.end())
		{
			//ignore pairs with bad estimates
			if(!pri->second->error && !post->second->error)
			{
				//dl/(dV_A)^2
				priNum = pri->second->f2j + 0.5*pri->second->f1j;
				postNum = post->second->f2j + 0.5*post->second->f1j;
				priDenom = 1 + priNum*varA + pri->second->f2j*varD;
				postDenom = 1 + postNum*varA + post->second->f2j*varD;
				priDenom = priDenom*priDenom;
				postDenom = postDenom*postDenom;

				priorBit = (priNum*priNum)/priDenom;
				postBit = (postNum*postNum)/postDenom;

				ans[0][0] += postBit - priorBit;

				//dl/(dV_A dV_D)				
				priorBit = (priNum*pri->second->f2j)/priDenom;
				postBit = (postNum*post->second->f2j)/postDenom;

				ans[0][1] += postBit - priorBit;

				//dl/(dV_D)^2				
				priorBit = (pri->second->f2j*pri->second->f2j)/priDenom;
				postBit = (post->second->f2j*post->second->f2j)/postDenom;

				ans[1][1] += postBit - priorBit;
			};
		};

		++pri;
		if(post != posteriors->snpIBDEstimates.end()) ++post;
	}while(pri != priors->snpIBDEstimates.end());

	//The Hessian matrix is symetric so only calculate half and then copy	
	ans[1][0] = ans[0][1];	
	
	//setup the matrix with calculated values
	map<unsigned int, double> aCol;	

	for(unsigned int col = 0; col < 2; ++col)
	{		
		for(unsigned int row = 0; row < 2; ++row)
		{
			aCol[row + 1] = ans[col][row];	
		};
		hessianMatrix[col + 1] = aCol;
	};

};

//! Returns derivative with respect to one parameter.
double ModelIBD::getDeriv(unsigned int & paraToFit) const
{
	double ans = 0;
	double priorBit, postBit;
	double priNum, postNum;

	//loop through prior and posterior estimates and calculate neglog like
	map<string, IBDData *>::const_iterator pri = priors->snpIBDEstimates.begin();
	map<string, IBDData *>::const_iterator post = posteriors->snpIBDEstimates.begin();

	do{
		//check prior and posteriors are for the same pair
		if(post == posteriors->snpIBDEstimates.end())
		{
			//find correct pair in posteriors
			post = posteriors->snpIBDEstimates.find(pri->first);
		}
		else if(post != posteriors->snpIBDEstimates.end() && pri->first != post->first)
		{
			//find correct pair in posteriors
			post = posteriors->snpIBDEstimates.find(pri->first);			
		};

		if(post != posteriors->snpIBDEstimates.end())
		{
			//ignore pairs with bad estimates
			if(!pri->second->error && !post->second->error)
			{
				priNum = pri->second->f2j + 0.5*pri->second->f1j;
				postNum = post->second->f2j + 0.5*post->second->f1j;

				if(paraToFit == 1)
				{
					//dl/dV_A					
					priorBit = priNum/(1 + priNum*varA + pri->second->f2j*varD);
					postBit = postNum/(1 + postNum*varA + post->second->f2j*varD);

					ans += priorBit - postBit;
				}
				else
				{
					//dl/dV_D				
					priorBit = pri->second->f2j/(1 + priNum*varA + pri->second->f2j*varD);
					postBit = post->second->f2j/(1 + postNum*varA + post->second->f2j*varD);

					ans += priorBit - postBit;
				};
			};
		};

		++pri;
		if(post != posteriors->snpIBDEstimates.end()) ++post;
	}while(pri != priors->snpIBDEstimates.end());

	return ans;
};

//! Returns 2nd derivative with respect to one parameter.
double ModelIBD::get2ndDeriv(unsigned int & paraToFit) const
{
	double ans = 0;
	
	double priorBit, postBit;
	double priNum, postNum, priDenom, postDenom;

	//loop through prior and posterior estimates and calculate neglog like
	map<string, IBDData *>::const_iterator pri = priors->snpIBDEstimates.begin();
	map<string, IBDData *>::const_iterator post = posteriors->snpIBDEstimates.begin();

	do{
		//check prior and posteriors are for the same pair
		if(post == posteriors->snpIBDEstimates.end())
		{
			//find correct pair in posteriors
			post = posteriors->snpIBDEstimates.find(pri->first);
		}
		else if(post != posteriors->snpIBDEstimates.end() && pri->first != post->first)
		{
			//find correct pair in posteriors
			post = posteriors->snpIBDEstimates.find(pri->first);			
		};

		if(post != posteriors->snpIBDEstimates.end())
		{
			//ignore pairs with bad estimates
			if(!pri->second->error && !post->second->error)
			{
				priNum = pri->second->f2j + 0.5*pri->second->f1j;
				postNum = post->second->f2j + 0.5*post->second->f1j;
				priDenom = 1 + priNum*varA + pri->second->f2j*varD;
				postDenom = 1 + postNum*varA + post->second->f2j*varD;
				priDenom = priDenom*priDenom;
				postDenom = postDenom*postDenom;

				if(paraToFit == 1)
				{
					//dl/(dV_A)^2
					priorBit = (priNum*priNum)/priDenom;
					postBit = (postNum*postNum)/postDenom;

					ans += postBit - priorBit;
				}
				else
				{		
					//dl/(dV_D)^2				
					priorBit = (pri->second->f2j*pri->second->f2j)/priDenom;
					postBit = (post->second->f2j*post->second->f2j)/postDenom;

					ans += postBit - priorBit;
				};
		
			};
		};

		++pri;
		if(post != posteriors->snpIBDEstimates.end()) ++post;
	}while(pri != priors->snpIBDEstimates.end());


	return ans;	
};
