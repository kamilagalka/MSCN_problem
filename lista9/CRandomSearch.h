#pragma once
#include "CMscnProblem.h"
#include "COptimizer.h"

class CRandomSearch:public COptimizer
{
public:
	CRandomSearch();
	CRandomSearch(CProblem* cProblem,int _numOfTries);
	~CRandomSearch() { delete problem; }

	Array<double>* search(int numOfTries, int& errCode);
	double* search();

private:
	int numOfTries;
	CProblem *problem;
};