#pragma once
#include "Array.h"
#include <iostream>

class CProblem {
public:
	virtual double dGetQuality(double* _solution, int& errCode) { 
		//std::cout << "\nTu niestety\n";
		return 10; };


	virtual Array<double>* generateSolution(int& errCode) { 
		//std::cout << "\nTu niestety\n\n\n";

		return new Array<double>; };


	virtual bool bConstraintsSatisfied(double* solution, int& errCode) { 
		//std::cout << "\nTu niestety\n\n\n";
		return false; };


	virtual int getSolutionSize() { 
	//	std::cout << "\nTu niestety\n\n\n";
		return 0; 
	};
};