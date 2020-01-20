#pragma once
#include "Array.h"
#include <iostream>

class CProblem {
public:
	virtual double dGetQuality(double* _solution, int& errCode) = 0;
	virtual Array<double>* generateSolution(int& errCode) = 0;
	virtual bool bConstraintsSatisfied(double* solution, int& errCode) = 0;
	virtual int getSolutionSize() = 0;
};