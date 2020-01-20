#pragma once
#include <iostream>
#include "Matrix.h"
#include "Array.h"
#include "CProblem.h"

class CMscnProblem : public CProblem{
public:
	CMscnProblem(int iD, int iF, int iM, int iS);
	~CMscnProblem();
	void vSetD(int iVal, int& successInfo);
	void vSetF(int iVal, int& successInfo);
	void vSetM(int iVal, int& successInfo);
	void vSetS(int iVal, int& successInfo);

	double dGetQuality(double *pdSolution, int& errCode);
	bool bConstraintsSatisfied(double *pdSolution, int& isCorrect);

	double getMinAt(double *pdSolution, int index,int& errCode);
	double getMaxAt(double *pdSolution, int index,int& errCode);

	void readProblemFile(std::string sFilename, int& errCode);
	void saveProblemFile(std::string sFilename, int& errCode);

	void readSolutionFile(std::string sFilename, int&errCode);
	void saveSolutionFile(std::string sFilename, int&errCode);

	void vGenerateInstance(int iInstanceSeed, int& errCode);
	Array<double>* generateSolution(int& errCode);
	Array<double>* repairSolution(double* pdSolution,int & errCode);

	Array<double>* getSolution(int& errCode);
	void setSolutionArray(Array<double>* newSolution,int& errCode);
	void setSolution(double* newSolution);

	int getSolutionSize();
	void printProblem();

private:
	int D;
	int F;
	int M;
	int S;

	Matrix<double> *cd;
	Matrix<double> *cf;
	Matrix<double> *cm;

	Array<double> *sd;
	Array<double> *sf;
	Array<double> *sm;
	Array<double> *ss;

	Matrix<double> *xd;
	Matrix<double> *xf;
	Matrix<double> *xm;

	Array<double> *ud;
	Array<double> *uf;
	Array<double> *um;

	Array<double> *ps;

	Matrix<double> *xdminmax;
	Matrix<double> *xfminmax;
	Matrix<double> *xmminmax;

	Array<double>* solution;
};