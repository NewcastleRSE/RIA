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


/*! \file main.h
    \brief This file defines a global variable to suppress output to screen. 
    
*/

// Uncomment if trying to compile using WINDOWS
#define USING_WINDOWS

#ifndef __MAIN
#define __MAIN

extern bool outputToScreen; 
extern ofstream logFile; 


template<typename T>
//! Outputs message to screen and log file
void out(const T & text)
{	
	if(outputToScreen) cout << text;
	logFile << text;
};

template<typename T>
//! Outputs to screen only.
void screenOut(const T & text)
{	
	if(outputToScreen) cout << text;	
};

template<typename T>
void screenFlush(const T & text)
{	
	if(outputToScreen) cout << flush;
};

template<typename T>
//! Outputs error message to screen and log file
void outErr(const T & text)
{
	cerr << text;
	logFile << text;
};

#endif
