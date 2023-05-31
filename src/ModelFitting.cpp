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


/*! \file ModelFitting.cpp
    \brief This file contains the source of the model fitting methods.
    
*/

#include <string>
#include <map>
#include <set>
#include <list>
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>

using namespace std; // initiates the "std" or "standard" namespace

#include "ModelFitting.h"
#include "Model.h"
#include "main.h"

//! Sets parameter values in map copied from another map.
void setParameterValues(map<unsigned int, double> & parametersToSet, map<unsigned int, double> & values)
{
	for(map<unsigned int, double>::const_iterator j = values.begin();  j != values.end(); ++j)
	{			
		parametersToSet[j->first] = j->second;
	};
};

//! Returns vector solution, ans, to matrix equation, X*ans = v, X is n by n and ans and v are n by 1.
map<unsigned int, double> getSolnMatrixEqun(map<unsigned int, map<unsigned int, double> > & matrix, map<unsigned int, double> & vect)
{
	map<unsigned int, double> ans = vect;

	//loop thro' rows of matrix, m, and the inverse, i
	map<unsigned int, map<unsigned int, double> >::iterator mrow = matrix.begin();	
	map<unsigned int, map<unsigned int, double> >::iterator mrow2;	
	map<unsigned int, double>::iterator mcol;	
	map<unsigned int, double>::iterator mcol2;	
	map<unsigned int, double>::iterator ansIt = ans.begin();
	map<unsigned int, double>::iterator ansIt2;
	
	double factor;
	unsigned int rowNo = 1;
	unsigned int colNo = 1;

	for( ; mrow != matrix.end(); ++mrow, ++ansIt)
	{
		//set column to the first column in the row
		mcol = mrow->second.begin();		
		colNo = 1;

		//advance the column until the the row no. is equal to the column no.
		while(colNo != rowNo)
		{
			mcol++;			
			++colNo;
		};

		//divide the row in (m and i) by the value in the matrix, m, at (rowNo, colNo)
		factor = 1.0/(mcol->second); 
		mcol->second = 1;		
		mcol++;		

		//scale the remaining elements in the row - if there are any
		while(mcol != mrow->second.end())
		{
			mcol->second *= factor;			
			mcol++;			
		};

		//scale answer vector by the same factor
		ansIt->second *= factor;

		//subtract the row in question away from the remaining rows scaled  s.t. column = mrow will be zero below this row in matrix m
		mrow2 = mrow;
		ansIt2 = ansIt;		
		mrow2++;
		ansIt2++;
		
		//loop thro' remaining rows
		while(mrow2 != matrix.end())
		{
			//set column iterators to begining of the rows
			mcol2 = mrow2->second.begin();			
			mcol = mrow->second.begin();			

			//advance column of matrix, m, to the same as the main row being considered, with value rowNo
			colNo = 1;
			while(colNo != rowNo)
			{
				mcol++;
				mcol2++;				
				++colNo;
			};

			factor = mcol2->second; //factor to multiple row, rowNo by to take away from subseq. rows
			mcol2->second -= (mcol->second)*factor;//0;
			mcol++;
			mcol2++;

			//subtract scaled row for the rest of the matrix, m
			while(mcol2 != mrow2->second.end())
			{
				mcol2->second -= (mcol->second)*factor;				
				mcol++;
				mcol2++;				
			};

			//subtract scaled value from ans vector
			ansIt2->second -= (ansIt->second)*factor;			

			mrow2++;
			ansIt2++;			
		};//end of performing row operations to set column (=rowNo) to zeroes below row = rowNo

		++rowNo;	
	};//end of performing row operations to set lower left of matrix, m, to zeroes


	//Now reduce the upper right of matrix, m, to zero
	map<unsigned int, map<unsigned int, double> >::reverse_iterator mrowre = matrix.rbegin();	
	map<unsigned int, map<unsigned int, double> >::reverse_iterator mrowre2 = matrix.rbegin();	
	map<unsigned int, double>::reverse_iterator ansItre = ans.rbegin();
	map<unsigned int, double>::reverse_iterator ansItre2;
	map<unsigned int, double>::reverse_iterator mcolre2;
	map<unsigned int, double>::reverse_iterator mcolre;

	rowNo = matrix.size();

	for( ; mrowre != matrix.rend(); ++mrowre, ++ansItre)
	{

		mrowre2 = mrowre;		
		ansItre2 = ansItre;
		mrowre2++;
		ansItre2++;		

		//loop tho' the remaining rows backwards - if there are any
		while(mrowre2 != matrix.rend())
		{			
			//set column iterators to begining of the rows
			mcolre2 = mrowre2->second.rbegin();			
			mcolre = mrowre->second.rbegin();		

			//advance column of matrix, m, to the same as the main row being considered, with value rowNo
			colNo = mrowre2->second.size();//size will be 4
			while(colNo != rowNo)
			{
				mcolre++;
				mcolre2++;				
				--colNo;
			};

			factor = mcolre2->second; //factor to multiple row, rowNo by to take away from subseq. rows
			mcolre2->second -= (mcolre->second)*factor;//0;
			mcolre++;
			mcolre2++;

			//subtract scaled row from the rest of the matrix, m
			while(mcolre2 != mrowre2->second.rend())//could stop at when col < row
			{
				mcolre2->second -= (mcolre->second)*factor;				
				mcolre++;
				mcolre2++;				
			};

			//subtract scaled value from ans vector
			ansItre2->second -= (ansItre->second)*factor;			

			mrowre2++;
			ansItre2++;			
		};

		--rowNo;
	};

	return ans;
};

//! Function used in Newton's Method, sum of square distance between 2 vectors.
double getDistanceVectors(map<unsigned int, double> & p1, map<unsigned int, double> & p2)
{
	double ans = 0;
	map<unsigned int, double>::const_iterator v1 = p1.begin();
	map<unsigned int, double>::const_iterator v2 = p2.begin();

	while(v1 != p1.end())	
	{
		ans += (v1->second - v2->second)*(v1->second - v2->second);
		++v1; ++v2;
	};

	return ans;
};

//! Finds the minimum negative log likelihood.
bool FindFit::newtonsMethod(double & eval, set<unsigned int> & parasToFit, double accuracy, unsigned int maxIterations) const
{
	double distance = 0;
	bool fitVarA = false;
	bool fitVarD = false;
	set<unsigned int>::const_iterator fb = parasToFit.find(1);
	if(fb != parasToFit.end()) fitVarA = true;
	fb = parasToFit.find(2);
	if(fb != parasToFit.end()) fitVarD = true;
	
	map<unsigned int, double> prevPoint;	
	map<unsigned int, double> nextPoint;

	unsigned int noIterations = 1;

	//setup prevPoint for the initial point	
	map<unsigned int, double> modelParameters = modelToFit->getParameters();

	for(map<unsigned int, double>::const_iterator mp = modelParameters.begin(); mp != modelParameters.end(); ++mp)
	{
		if(parasToFit.find(mp->first) != parasToFit.end())
		{
			prevPoint[mp->first] = mp->second;
			nextPoint[mp->first] = mp->second;
		};
	};
 
	do{
		getNextPoint(nextPoint, fitVarA, fitVarD);

		modelToFit->setNewParameters(nextPoint);

		//check for convergence compare last two points and/or loglikelihoods
		distance = getDistanceVectors(prevPoint, nextPoint);

		if(distance*0 != 0) return false;
		else if(distance < accuracy) break;

		setParameterValues(prevPoint, nextPoint);
		noIterations++;
	}while(noIterations <= maxIterations);

	eval = evaluate();

	return (noIterations <= maxIterations) && (eval*0 == 0) && (eval <= 0);
};

//! Calculates the next point for one iteration of Newton's method for optimization.
void FindFit::getNextPoint(map<unsigned int, double> & point, const bool & fitVarA, const bool & fitVarD) const
{
	map<unsigned int, double> pointChange = getPointChange(fitVarA, fitVarD);
	map<unsigned int, double>::const_iterator pc = pointChange.begin();
	map<unsigned int, double>::const_iterator pp = point.begin();

	for( ; pp != point.end(); ++pp, ++pc)
	{
		point[pp->first] -= pc->second;
	};

};

//! Returns vector on how to change the best fit parameters when using Newton's Method.
map<unsigned int, double> FindFit::getPointChange(const bool & fitVarA, const bool & fitVarD) const
{	
	map<unsigned int, double> ans;

	map<unsigned int, map<unsigned int, double> > hessianMatrix; // col, row
	map<unsigned int, double> gradientVector;
	map<unsigned int, double> aVector;

	//get the gradient vector, 1st derivatives with repect to each parameters
	modelToFit->getGradientVector(gradientVector, fitVarA, fitVarD);

	//get the Hessian matrix, 2nd derivatives with repect to each pair of parameters
	modelToFit->getHessianMatrix(hessianMatrix, fitVarA, fitVarD);

	//solve Hess*ans = gradVect
	ans = getSolnMatrixEqun(hessianMatrix, gradientVector);
	return ans;
};

//! Fits Newton's for one variable.
bool FindFit::newtonsMethodForOneVariable(double & eval, unsigned int & paraToFit, double accuracy, unsigned int maxSteps) const
{
	map<unsigned int, double> parametersFit = modelToFit->getParameters();
	double initialValue = evaluate();
	double initialParaValue = modelToFit->getParameter(paraToFit);
	double parameterValue = initialParaValue;
	double paraStep;
	unsigned int steps = 0;

	while(true)
	{
		paraStep = modelToFit->getDeriv(paraToFit) / modelToFit->get2ndDeriv(paraToFit);
		parameterValue = parameterValue - paraStep;
		
		modelToFit->setParameter(paraToFit, parameterValue);

		if(++steps >= maxSteps)	break;
		else if(paraStep*0 != 0 || fabs(paraStep) < accuracy) break;
	};

	eval = evaluate();
	if(eval > initialValue || eval*0 != 0 || eval > 0)
	{		
		return false;
	};

	return true;
};

//! Returns the sum of parameters.
map<unsigned int, double> getSumParameters(map<unsigned int, map<unsigned int, double> > & simplex)
{
	map<unsigned int, double> ans;
	for(map<unsigned int, double>::const_iterator p = simplex.begin()->second.begin(); p != simplex.begin()->second.end(); ++p)
	{	
			ans[p->first] = 0;
	};

	for(map<unsigned int, map<unsigned int, double> >::const_iterator vertex = simplex.begin(); vertex != simplex.end(); ++vertex)
	{		
		for(map<unsigned int, double>::const_iterator para = vertex->second.begin(); para != vertex->second.end(); ++para)
		{			
			ans[para->first] += para->second;
		};		 
	};

	return ans;
};

//! Fits model using the downhill simplex method.
bool FindFit::downhillSimplex(double & eval, set<unsigned int> parasToFit, double accuracy, unsigned int maxIterations,  double e3) const
{
	bool noStoppingOn3 = false;
	if(accuracy == 0) {noStoppingOn3 = true; accuracy = 1e-6;};

	double aNo;
	map<unsigned int, double> parametersFit, parametersAve, parametersRef, parametersNew, parametersExpan;
	map<unsigned int, map<unsigned int, double> > simplex; //simplex of points, 0 is origin (initial guess), plus N (no. of parameters) points

	map<unsigned int, double> modelParameters = modelToFit->getParameters();
	set<unsigned int> parametersToFit;
	if(parasToFit.empty())
	{		
		for(map<unsigned int, double>::const_iterator pf = modelParameters.begin(); pf != modelParameters.end(); ++pf) {parametersToFit.insert(pf->first);};
	}
	else
	{		
		parametersToFit = parasToFit;
	};

	//setup simplex
	map<unsigned int, double> initGuessParas;
	
	for(map<unsigned int, double>::const_iterator mp = modelParameters.begin(); mp != modelParameters.end(); ++mp)
	{
		if(parametersToFit.find(mp->first) != parametersToFit.end())
		{
			initGuessParas[mp->first] = mp->second;
		};
	};

	map<unsigned int, double> vertexEvaluations;
	simplex[0] = initGuessParas;
	modelToFit->setNewParameters(initGuessParas);
	vertexEvaluations[0] = evaluate();

	unsigned int highest;
	unsigned int secondHighest;
	unsigned int lowest;
	unsigned int stoppingReason;
	double accuracy2 = accuracy*accuracy;
	double evalMaxOld, evalMax, eval2ndMax, evalRef, sum, evalNew, evalExpan, evalMin;	
	double num;
	double factor = 1.1;

	for(map<unsigned int, double>::const_iterator igp = initGuessParas.begin(); igp != initGuessParas.end(); ++igp)
	{
		setParameterValues(simplex[igp->first], initGuessParas);
		num = factor*(igp->second); //increase a bit, proportional to the size of the parameter		
		if(num == 0) num = 0.1;

		simplex[igp->first][igp->first] += num;

		//record evaluations at each vertex of the simplex
		modelToFit->setNewParameters(simplex.find(igp->first)->second);
		vertexEvaluations[igp->first] = evaluate();
	};

	double noParameters = (double)(simplex.size() - 1);	
	map<unsigned int, double> sumParameters;
	map<unsigned int, double>::const_iterator m,m2,av2,m3,r3,low,hi;

	unsigned int noIterations = 0;
	evalMaxOld = 0;

	for (;;)
	{
		//find the highest and the lowest vertex evaluation
		highest = 0; lowest = 0; secondHighest = 0; 
		for(map<unsigned int, double>::const_iterator ve = vertexEvaluations.begin(); ve != vertexEvaluations.end(); ++ve)
		{
			if(ve->second > vertexEvaluations.find(highest)->second) {secondHighest = highest; highest = ve->first;}
			else if(ve->second > vertexEvaluations.find(secondHighest)->second) secondHighest = ve->first;
			else if(ve->second < vertexEvaluations.find(lowest)->second) lowest = ve->first;
		};

		sumParameters = getSumParameters(simplex);

		//set average of parameters (except highest pt)
		m = simplex.find(highest)->second.begin();
		for(map<unsigned int, double>::const_iterator sp = sumParameters.begin(); sp != sumParameters.end(); ++sp, ++m)
		{			
			parametersAve[sp->first] = ((sp->second) - (m->second)) / noParameters;
		};

		//set reflection of simplex through the highest pt
		m2 = simplex.find(highest)->second.begin();
		for(map<unsigned int, double>::const_iterator av = parametersAve.begin(); av != parametersAve.end(); ++av, ++m2)
		{	
			aNo = (2 * (av->second)) - (m2->second);
			if(aNo < 0) parametersRef[av->first] = 0; //can only be negative if permitted
			else parametersRef[av->first] = aNo;
		};

		evalMax = vertexEvaluations.find(highest)->second;
		evalMin = vertexEvaluations.find(lowest)->second;
		eval2ndMax = vertexEvaluations.find(secondHighest)->second;
		modelToFit->setNewParameters(parametersRef); 
		evalRef = evaluate();

		if(evalMin > evalRef)
		{
			//Attempt Expansion
			av2 = parametersAve.begin();
			for(map<unsigned int, double>::const_iterator r = parametersRef.begin(); r != parametersRef.end(); ++r, ++av2)
			{				
				aNo = (2 * (r->second)) - (av2->second);
				if(aNo < 0) parametersExpan[r->first] = 0;
				else parametersExpan[r->first] = aNo;
			};
			modelToFit->setNewParameters(parametersExpan);
			evalExpan = evaluate();
			if(evalExpan <= evalRef) {evalNew = evalExpan; setParameterValues(parametersNew, parametersExpan);}
			else {evalNew = evalRef; setParameterValues(parametersNew, parametersRef);};
		}
		else if(eval2ndMax > evalRef)
		{
			//Do Reflection
			setParameterValues(parametersNew, parametersRef); evalNew = evalRef;
		}
		else
		{
			//Perform Contraction
			if(evalMax <= evalRef)
			{
				m3 = simplex.find(highest)->second.begin();
				for(map<unsigned int, double>::const_iterator av3 = parametersAve.begin(); av3 != parametersAve.end(); ++av3, ++m3)
				{			
					parametersNew[av3->first] = ((m3->second) + (av3->second)) / 2 ;
				};
			}
			else
			{
				r3 = parametersRef.begin();
				for(map<unsigned int, double>::const_iterator av4 = parametersAve.begin(); av4 != parametersAve.end(); ++av4, ++r3)
				{			
					parametersNew[av4->first] = ((r3->second) + (av4->second)) / 2 ;
				};

			};

			modelToFit->setNewParameters(parametersNew); 
			evalNew = evaluate();
		};

		if(evalNew > evalMax)
		{
			//Shrink the simplex toward the lowest vertex
			for(map<unsigned int, map<unsigned int, double> >::iterator sim = simplex.begin(); sim != simplex.end(); ++sim)
			{
				if(sim->first != lowest)
				{
					low = simplex.find(lowest)->second.begin();
					for(map<unsigned int, double>::iterator vert = sim->second.begin(); vert != sim->second.end(); ++vert, ++low)
					{						
						vert->second = ((vert->second) + (low->second)) / 2;						
					};
					modelToFit->setNewParameters(sim->second);
					vertexEvaluations[sim->first] = evaluate();
				};
			};
		};

		//test for stopping criteria, if new vertex has not changed much then stop		
		sum = 0;
		hi = simplex.find(highest)->second.begin();
		for(map<unsigned int, double>::const_iterator pn = parametersNew.begin(); pn != parametersNew.end(); ++pn, ++hi)
		{			

			sum +=  ((pn->second) - (hi->second)) * ((pn->second) - (hi->second));
		};

		if(sum < accuracy2)
		{
			stoppingReason = 1; break;
		}
		else if(fabs(evalMaxOld - evalMax) < e3 && !noStoppingOn3 && noIterations > 1)
		{
			stoppingReason = 3; break;
		}
		else if(noIterations++ >= maxIterations)
		{
			stoppingReason = 2; break;
		};
		
		evalMaxOld = evalMax;

		//form the new vertex by replacing the maximum vertex with the new vertex
		setParameterValues(simplex[highest], parametersNew);
		vertexEvaluations[highest] = evalNew;

	};//end of main loop

	modelToFit->setNewParameters(simplex.find(lowest)->second);
	
	setParameterValues(parametersFit, simplex.find(lowest)->second);

	eval = evalMin;

	return (eval*0 == 0) && (stoppingReason != 2) && (eval <= 0);
};
