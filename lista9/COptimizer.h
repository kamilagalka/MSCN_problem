#pragma once
#include "CProblem.h"

class COptimizer
{
public:
	virtual double* search() = 0;
private:
	CProblem* problem;
};