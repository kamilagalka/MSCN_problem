// lista9.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "pch.h"
#include "CMscnProblem.h"
#include "ErrorCodes.h"
#include "CRandom.h"
#include "CMscnProblem.h"
#include "CRandomSearch.h"
#include "CDiffEvol.h"
#include "CTimer.h"
#include <iostream>


//-----------------MAIN DO ZAD 3---------------
int main() {
	CTimer timer;
	timer.start();
	
	int errCode = IS_OK;
	CMscnProblem* problem;
	problem = new CMscnProblem(1, 1, 1, 1);
	problem->readProblemFile("plikProblemu2.txt", errCode);

	CDiffEvol diffEvol(problem, 100,5);
	double* solutionDiffEvol = diffEvol.search();
	std::cout << "\nSOLUTION: \n";
	for (int i = 0; i < problem->getSolutionSize(); i++) {
		std::cout << "  " << solutionDiffEvol[i];
	}
	std::cout << "\nQUALITY: \n";
	std::cout << problem->dGetQuality(solutionDiffEvol, errCode);
	std::cout << "\nCONSTRAINTS SATISFIED: \n";
	std::cout << problem->bConstraintsSatisfied(solutionDiffEvol, errCode);


	problem->setSolution(solutionDiffEvol);


	timer.stop();


	std::cout << "\n\n\nTime: ";
	std::cout << timer.getSecs();
}

//--------------MAIN DO GENEROWANIA ROZWIAZAN------------------------
/*
int main() {
	int errCode = IS_OK;
	CMscnProblem* problem;
	problem = new CMscnProblem(1, 1, 1, 1);
	problem->readProblemFile("plikProblemu2.txt", errCode);

	CRandomSearch rsearch(problem, 1000);
	//problem->setSolutionArray(rsearch.search(10000,errCode),errCode);
	problem->setSolution(rsearch.search());

	std::cout << "\nSOLUTION: \n";
	problem->getSolution(errCode)->printArray();
	std::cout << "\nQUALITY: \n";
	std::cout << problem->dGetQuality(problem->getSolution(errCode)->getTable(), errCode);
	std::cout << "\nCONSTRAINTS SATISFIED: \n";
	std::cout << problem->bConstraintsSatisfied(problem->getSolution(errCode)->getTable(), errCode);

	std::cout << "\n\n-------------------------------------------------------------\n\n";

	CDiffEvol diffEvol(problem, 100);
	double* solutionDiffEvol = diffEvol.search();
	std::cout << "\nSOLUTION: \n";
	for (int i = 0; i < problem->getSolutionSize(); i++) {
		std::cout << "  " << solutionDiffEvol[i];
	}
	std::cout << "\nQUALITY: \n";
	std::cout << problem->dGetQuality(solutionDiffEvol, errCode);
	std::cout << "\nCONSTRAINTS SATISFIED: \n";
	std::cout << problem->bConstraintsSatisfied(solutionDiffEvol, errCode);


	problem->setSolution(solutionDiffEvol);
	problem->saveSolutionFile("solution.txt", errCode);

}
*/


//-------------------MAIN DO WERYFIKACJI ROZWIAZAN-------------------------
/*

int main() {
	int errCode = IS_OK;
	CMscnProblem* problem;
	problem = new CMscnProblem(1, 1, 1, 1);
	problem->readProblemFile("plikProblemu2.txt", errCode);
	problem->readSolutionFile("plikRozwiazania2.txt", errCode);

	std::cout << "\n\nQuality:  ";
	std::cout<<problem->dGetQuality(problem->getSolution(errCode)->getTable(), errCode);

	std::cout << "\n\nConstraints: ";
	std::cout << problem->bConstraintsSatisfied(problem->getSolution(errCode)->getTable(), errCode);

	std::cout << "\n\n\n";
	problem->getSolution(errCode)->printArray();
}
*/



/*
int main() {
	int errCode = IS_OK;
	CMscnProblem* problem;
	problem = new CMscnProblem(2,2,2,2);
	problem->vGenerateInstance(5, errCode);
	problem->printProblem();
	std::cout << "\n\n-------------------------------------------------------------\n\n";
	
	CRandomSearch rsearch(problem,1000);
	//problem->setSolutionArray(rsearch.search(10000,errCode),errCode);
	problem->setSolution(rsearch.search());
	
	std::cout << "\nSOLUTION: \n";
	problem->getSolution(errCode)->printArray();
	std::cout << "\nQUALITY: \n";
	std::cout << problem->dGetQuality(problem->getSolution(errCode)->getTable(), errCode);
	std::cout << "\nCONSTRAINTS SATISFIED: \n";
	std::cout<<problem->bConstraintsSatisfied(problem->getSolution(errCode)->getTable(),errCode);

	std::cout << "\n\n-------------------------------------------------------------\n\n";
	
	CDiffEvol diffEvol(problem, 50);
	double* solutionDiffEvol = diffEvol.search();
	std::cout << "\nSOLUTION: \n";
	for (int i = 0; i < problem->getSolutionSize(); i++) {
		std::cout << "  " << solutionDiffEvol[i];
	}
	std::cout << "\nQUALITY: \n";
	std::cout << problem->dGetQuality(solutionDiffEvol, errCode);
	std::cout << "\nCONSTRAINTS SATISFIED: \n";
	std::cout << problem->bConstraintsSatisfied(solutionDiffEvol, errCode);


}*/