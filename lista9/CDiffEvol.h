#pragma once
#include "CMscnProblem.h"
#include "Array.h"
#include "Matrix.h"
#include "CRandom.h"
#include "CProblem.h"
#include "COptimizer.h"
#include "CTimer.h"

class CDiffEvol:COptimizer {
public:
	CDiffEvol(CMscnProblem* _problem,int populationSize,double time);
	~CDiffEvol();
	void initPopulation();
	bool checkStopCondition(int counter);
	bool checkStopCondition(double timeInSecs);
	double* getRandomInd();
	bool individualsAreDifferent(double* ind, double* baseInd, double* addInd0, double* addInd1);
	double* search();

	void saveQualitiesToFile();

private:
	CTimer timer;
	double time;
	double** population;
	int populationSize;
	int genotypeSize;
	CProblem* problem;
};