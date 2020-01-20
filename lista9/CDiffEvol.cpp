#include "pch.h"
#include "CDiffEvol.h"
#include <iostream>
#include <fstream>



CDiffEvol::CDiffEvol(CMscnProblem * _problem, int _populationSize,double _time)
{
	problem = _problem;
	populationSize = _populationSize;
	genotypeSize = problem->getSolutionSize();
	time = _time;
}

CDiffEvol::~CDiffEvol()
{
}

void CDiffEvol::initPopulation()
{
	int errCode=0;
	population = new double*[populationSize];
	double best = 0;
	for (int i = 0; i < populationSize; i++) {
		population[i] = (problem->generateSolution(errCode))->getTable();
	}


}

bool CDiffEvol::checkStopCondition(int counter)
{
	return counter == 1000;
}

bool CDiffEvol::checkStopCondition(double timeInSecs)
{
	return timeInSecs >= time;
}

double * CDiffEvol::getRandomInd()
{
	CRandom random;
	int index = random.getRandomInt(0,populationSize-1);
	return population[index];
}

bool CDiffEvol::individualsAreDifferent(double * ind, double * baseInd, double * addInd0, double * addInd1)
{
	for (int i = 0; i < problem->getSolutionSize(); i++) {
		if (ind[i] != baseInd[i] || ind[i] != addInd0[i] || addInd0[i] != addInd1[i]) {
			return true;
		}
		if (baseInd[i] != addInd0[i] || baseInd[i] != addInd1[i]) {
			return true;
		}
		if (addInd0[i] != addInd1[i]) {
			return true;
		}

	}
	return false;
}

double* CDiffEvol::search()
{
	initPopulation();
	int max = 0;
	int counter = 0;
	int errCode = 0;

	double diffWeight = 0.2;
	double crossProb = 0.8;

	CRandom random;
	double* ind;
	double* baseInd;
	double* addInd0;
	double* addInd1;

	double* result = new double[genotypeSize];

	double* indNew;
	indNew = new double[genotypeSize];
	timer.start();

	while (checkStopCondition(timer.getSecs()) == false) {
		for (int i = 0; i < populationSize; i++) {
			ind = population[i];
			baseInd = getRandomInd();
			addInd0 = getRandomInd();
			addInd1 = getRandomInd();
			if (individualsAreDifferent(ind, baseInd, addInd0, addInd1)) {
				indNew = new double[genotypeSize];
				for (int geneOffset = 0; geneOffset < genotypeSize; geneOffset++) {
					if (random.getRandomDouble(0, 1) < crossProb) {
						indNew[geneOffset] = baseInd[geneOffset] + diffWeight * (addInd0[geneOffset] - addInd1[geneOffset]);
					}
					else {
						indNew[geneOffset] = ind[geneOffset];
					}
				}

				if (problem->dGetQuality(indNew, errCode) >= problem->dGetQuality(ind, errCode)&& problem->bConstraintsSatisfied(indNew,errCode)&&population[i]!=indNew) {
					population[i] = indNew;
		
					if (problem->dGetQuality(indNew, errCode) > max &&problem->bConstraintsSatisfied(indNew, errCode)) {
						max = problem->dGetQuality(indNew, errCode);
						errCode = 0;
						for (int l = 0; l < genotypeSize; l++) {
							result[l] = indNew[l];
						}
					}
				}

			}
		}
		counter++;
		timer.stop();
		//saveQualitiesToFile();
	}

	
	return result;
}
void CDiffEvol::saveQualitiesToFile()
{

	int errCode = 0;
	std::ofstream file;
	file.open("differentialEvolution.txt", std::ios::app);
	if (file.is_open())
	{
		for (int i = 0; i < populationSize - 1; i++)
		{
			file << problem->dGetQuality(population[i], errCode);
			file << ',';
		}
		file << problem->dGetQuality(population[populationSize - 1], errCode);
		file << "\n";
	}
}