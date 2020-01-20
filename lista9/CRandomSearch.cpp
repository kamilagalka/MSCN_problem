#include "pch.h"
#include "CRandomSearch.h"

CRandomSearch::CRandomSearch()
{
	problem = new CMscnProblem(1, 1, 1, 1);
	numOfTries = 100;
}

CRandomSearch::CRandomSearch(CMscnProblem * cProblem,int _numOfTries)
{
	problem = cProblem;
	numOfTries = _numOfTries;

}

Array<double>* CRandomSearch::search(int hereNumOfTries, int& errCode)
{
	if (hereNumOfTries <= 0) {
		errCode = INCORRECT_VALUE;
	}
	int improvements = 0;
	Array<double> *bestSolution = problem->generateSolution(errCode);
	double bestQuality = problem->dGetQuality(bestSolution->getTable(), errCode);
	bool bestConstraintsSatisfied = problem->bConstraintsSatisfied(bestSolution->getTable(), errCode);

	Array<double>* currentSolution;
	double currentQuality;
	bool currentConstraintsSatisfied = false;
	for (int i = 0; i < hereNumOfTries; i++) {
		currentSolution = problem->generateSolution(errCode);
		currentQuality = problem->dGetQuality(currentSolution->getTable(), errCode);
		currentConstraintsSatisfied = problem->bConstraintsSatisfied(currentSolution->getTable(), errCode);

		if (currentQuality > bestQuality || currentConstraintsSatisfied == true && bestConstraintsSatisfied == false) {
			improvements++;
			if (problem->bConstraintsSatisfied(currentSolution->getTable(), errCode)) {
				bestQuality = currentQuality;
				bestSolution = currentSolution;
			}
		}
	}
	
	std::cout << "\n\nNumber of improvements: " << improvements;
	std::cout << "\n\n";
	
	
	return bestSolution;
}
double* CRandomSearch::search()
{
	int errCode=0;
	return search(numOfTries, errCode)->getTable();
}