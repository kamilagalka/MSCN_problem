#pragma once
#include "CProblem.h"

class COptimizer
{
public:
	virtual double* search() { return 0; };
private:
	CProblem* problem;
};