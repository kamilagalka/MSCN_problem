#include "pch.h"
#include "CMscnProblem.h"
#include "ErrorCodes.h"
#include "Bounds.h"
#include "Matrix.h"
#include "Array.h"
#include "CRandom.h"

CMscnProblem::CMscnProblem(int iD, int iF, int iM, int iS)
{
	D = iD;
	F = iF;
	M = iM;
	S = iS;

	cd = new Matrix<double>(D, F);
	cf = new Matrix<double>(F, M);
	cm = new Matrix<double>(M, S);

	sd = new Array<double>(D);
	sf = new Array<double>(F);
	sm = new Array<double>(M);
	ss = new Array<double>(S);

	xd = new Matrix<double>(D, F);
	xf = new Matrix<double>(F, M);
	xm = new Matrix<double>(M, S);

	ud = new Array<double>(D);
	uf = new Array<double>(F);
	um = new Array<double>(M);

	ps = new Array<double>(S);


	xdminmax = new Matrix<double>(D, 2 * F);
	xfminmax = new Matrix<double>(F, 2 * M);
	xmminmax = new Matrix<double>(M, 2 * S);

	solution = new Array<double>(D*F + F * M + M * S);
}


CMscnProblem::~CMscnProblem()
{
	
	delete cd;
	delete cf;
	delete cm;

	delete sd;
	delete sf;
	delete sm;
	delete ss;


	delete ud;
	delete uf;
	delete um;

	delete xd;
	delete xf;
	delete xm;

	delete ps;

	delete xdminmax;
	delete xfminmax;
	delete xmminmax;
	
	delete solution;
	
}



void CMscnProblem::vSetD(int newD, int& errCode)
{
	if (newD <= 0) {
		errCode = NON_POSITIVE_VAL;
	}
	else if (newD != D) {
		cd->changeSizeX(newD, errCode);
		xd->changeSizeX(newD, errCode);
		sd->changeSize(newD, errCode);
		ud->changeSize(newD, errCode);
		xdminmax->changeSizeX(newD, errCode);

		D = newD;
		solution->changeSize(D*F + F * M + M * S,errCode);

	}
}

void CMscnProblem::vSetF(int newF, int& errCode)
{
	if (newF <= 0) {
		errCode = NON_POSITIVE_VAL;
	}
	else if (newF != F) {
		cd->changeSizeY(newF, errCode);
		cf->changeSizeX(newF, errCode);
		xd->changeSizeY(newF, errCode);
		xf->changeSizeX(newF, errCode);
		sf->changeSize(newF, errCode);
		uf->changeSize(newF, errCode);
		xdminmax->changeSizeY(2 * newF, errCode);
		xfminmax->changeSizeX(newF, errCode);

		F = newF;
		solution->changeSize(D*F + F * M + M * S, errCode);

	}
}

void CMscnProblem::vSetM(int newM, int& errCode)
{
	if (newM <= 0) {
		errCode = NON_POSITIVE_VAL;
	}
	else if (newM != M) {
		cf->changeSizeY(newM, errCode);
		cm->changeSizeX(newM, errCode);
		xf->changeSizeY(newM, errCode);
		xm->changeSizeX(newM, errCode);
		sm->changeSize(newM, errCode);
		um->changeSize(newM, errCode);
		xfminmax->changeSizeY(2 * newM, errCode);
		xmminmax->changeSizeX(newM, errCode);

		M = newM;
		solution->changeSize(D*F + F * M + M * S, errCode);

	}
}

void CMscnProblem::vSetS(int newS, int& errCode)
{
	if (newS <= 0) {
		errCode = NON_POSITIVE_VAL;
	}
	else if (newS != S) {
		cm->changeSizeY(newS, errCode);
		xm->changeSizeY(newS, errCode);
		ss->changeSize(newS, errCode);
		ps->changeSize(newS, errCode);
		xmminmax->changeSizeY(2 * newS, errCode);

		S = newS;
		solution->changeSize(D*F + F * M + M * S, errCode);

	}
}



double CMscnProblem::dGetQuality(double * pdSolution, int &errCode)
{
	int correctSize = D * F + F * M + M * S;
	double z = 0;

	if (pdSolution == NULL) {
		errCode = POINTER_IS_NULL;
	}

	for (int i = 0; i < correctSize; i++) {
		if (pdSolution[i] < 0) {
			errCode = NEGATIVE_VALUES;
			//std::cout << "WYWALAM BO NEGATIVE";
		}
	}
	errCode = 0;
	if (errCode == 0) {

		// ----- P -----
		double p = 0;
		int position = D * F + F * M;

		for (int i = 0; i < M; i++) {
			for (int j = 0; j < S; j++) {
				p += (ps->getValueAt(j, errCode))* pdSolution[position];
				position++;
			}
		}

		// ----- KU -----
		double ku = 0;
		position = 0;

		double sum = 0;
		for (int i = 0; i < D; i++) {
			for (int j = 0; j < F; j++) {
				sum += pdSolution[position];
				position++;
			}
			if (sum > 0) {
				ku += ud->getValueAt(i, errCode);
			}
			sum = 0;

		}
		sum = 0;
		for (int i = 0; i < F; i++) {
			for (int j = 0; j < M; j++) {

				sum += pdSolution[position];
				position++;
			}
			if (sum > 0) {

				ku += uf->getValueAt(i, errCode);
			}
			sum = 0;
		}

		sum = 0;
		for (int i = 0; i < M; i++) {
			for (int j = 0; j < S; j++) {
				sum += pdSolution[position];
				position++;
			}
			if (sum > 0) {

				ku += um->getValueAt(i, errCode);
			}
			sum = 0;

		}


		// ----- KT -----
		double kt = 0;
		position = 0;
		for (int i = 0; i < D; i++) {
			for (int j = 0; j < F; j++) {
				kt += (cd->getValueAt(i, j, errCode)) * pdSolution[position];
				position++;
			}
		}
		for (int i = 0; i < F; i++) {
			for (int j = 0; j < M; j++) {

				kt += (cf->getValueAt(i, j, errCode)) * pdSolution[position];
				position++;
			}
		}

		for (int i = 0; i < M; i++) {
			for (int j = 0; j < S; j++) {
				kt += (cm->getValueAt(i, j, errCode)) * pdSolution[position];
				position++;
			}
		}

		z = p - kt - ku;
	}
	errCode = 0;
	return z;

}

bool CMscnProblem::bConstraintsSatisfied(double * pdSolution, int& errCode)
{
	if (pdSolution == NULL) {
		errCode = POINTER_IS_NULL;
	}

	for (int i = 0; i < D*F + F * M + M * S; i++) {
		if (pdSolution[i] < 0) {
			errCode = NEGATIVE_VALUES;
		}
	}

	if (errCode == IS_OK) {
		double res = 0;
		int position = 0;


		for (int i = 0; i < D; i++) {
			for (int j = 0; j < F; j++) {
				res += pdSolution[position];
				position++;
			}
			if (res > sd->getValueAt(i, errCode)) {
				//std::cout << "\nwywalam w 1";
				return false;
			}
			res = 0;
		}


		for (int i = 0; i < F; i++) {
			for (int j = 0; j < M; j++) {
				res += pdSolution[position];
				position++;
			}
			if (res > sf->getValueAt(i, errCode)) {
				//std::cout << "\nwywalam w 2";

				return false;
			}
			res = 0;
		}


		for (int i = 0; i < M; i++) {
			for (int j = 0; j < S; j++) {
				res += pdSolution[position];
				position++;
			}

			if (res > sm->getValueAt(i, errCode)) {
				//std::cout << "\nwywalam w 3";

				return false;
			}
			res = 0;
		}

		position = D * F + F * M;
		for (int i = 0; i < S; i++) {
			position = D * F + F * M + i;
			for (int j = 0; j < M; j++) {
				res += pdSolution[position];
				position += S;
			}

			if (res > ss->getValueAt(i, errCode)) {
				//std::cout << "\nwywalam w 4";

				return false;

			}
			res = 0;
		}

		double sum1 = 0;
		int positionXf = D * F;

		for (int i = 0; i < F; i++) {
			for (int j = i; j < D*F; j += F) {
				sum1 += pdSolution[j];
			}
		}

		position = positionXf;
		double sum2 = 0;
		for (int i = 0; i < F; i++) {
			for (int j = 0; j < M; j++) {
				sum2 += pdSolution[position];
				position++;
			}
		}

		if (sum2 > sum1) {
			//std::cout << "\nwywalam w 5";

			return false;

		}

		sum1 = 0;
		sum2 = 0;

		int positionXm = D * F + F * M;

		for (int i = 0; i < M; i++) {
			for (int j = i; j < F*M; j += M) {
				sum1 += pdSolution[positionXf + j];
			}
		}

		position = positionXm;
		for (int i = 0; i < M; i++) {
			for (int j = 0; j < S; j++) {
				sum2 += pdSolution[position];
				position++;
			}
		}

		if (sum2 > sum1) {
			//std::cout << "\nsum1" << sum1;
			//std::cout << "\nsum" << sum2;
		//	std::cout << "\nwywalam w 6";

			return false;

		}
		return true;
	}
	else {
		//std::cout << "\nwywalam w 7";

		return false;
	}
}

double CMscnProblem::getMinAt(double * pdSolution, int index, int& errCode)
{
	int numOfAllElements = (D * F + F * M + M * S) - 1;
	int ims = (D * F + F * M) - 1;
	int ifm = (D * F) - 1;
	int idf = 0;
	if (index > ims) {
		return xmminmax->getValueAt(((index - (ims + 1)) / S), 2 * ((index - (ims + 1)) % S), errCode);
	}
	else if (index > ifm) {
		return xfminmax->getValueAt(((index - (ifm + 1)) / M), 2 * ((index - (ifm + 1)) % M), errCode);
	}
	else {
		return xdminmax->getValueAt(((index) / F), 2 * ((index%F)), errCode);
	}
	/*
	if (index >= D * F + F * M) {
		int position = 0;
		index = index - (D*F + F * M);
		for (int i = 0; i < M; i++) {
			for (int j = 0; j < 2 * S; j += 2) {
				if (position == index) {
					return xmminmax[i][j];
				}
				position++;

			}
		}
	}

	if (index >= D * F) {
		int position = 0;
		index = index - (D*F);
		for (int i = 0; i < F; i++) {
			for (int j = 0; j < 2 * M; j += 2) {
				if (position == index) {
					return xfminmax[i][j];
				}
				position++;

			}
		}
	}

	else {
		int position = 0;
		for (int i = 0; i < D; i++) {
			for (int j = 0; j < 2 * F; j += 2) {
				if (position == index) {
					return xdminmax[i][j];
				}
				position++;

			}
		}
	}
	*/



}


double CMscnProblem::getMaxAt(double * pdSolution, int index, int&errCode)
{
	int numOfAllElements = (D * F + F * M + M * S) - 1;
	int ims = (D * F + F * M) - 1;
	int ifm = (D * F) - 1;
	int idf = 0;
	if (index > ims) {
		return xmminmax->getValueAt(((index - (ims + 1)) / S), 2 * ((index - (ims + 1)) % S) + 1, errCode);
	}
	else if (index > ifm) {
		return xfminmax->getValueAt(((index - (ifm + 1)) / M), 2 * ((index - (ifm + 1)) % M) + 1, errCode);
	}
	else {
		return xdminmax->getValueAt(((index) / F), 2 * ((index%F)) + 1, errCode);
	}
	/*
	if (index >= D * F + F * M) {
		int position = 0;
		index = index - (D*F + F * M);
		for (int i = 0; i < M; i++) {
			for (int j = 1; j < 2 * S; j += 2) {
				if (position == index) {
					return xmminmax[i][j];
				}
				position++;

			}
		}
	}

	if (index >= D * F) {
		int position = 0;
		index = index - (D*F);
		for (int i = 0; i < F; i++) {
			for (int j = 1; j < 2 * M; j += 2) {
				if (position == index) {
					return xfminmax[i][j];
				}
				position++;

			}
		}
	}

	else {
		int position = 0;
		for (int i = 0; i < D; i++) {
			for (int j = 1; j < 2 * F; j += 2) {
				if (position == index) {
					return xdminmax[i][j];
				}
				position++;

			}
		}
	}
	*/
}


void CMscnProblem::readProblemFile(std::string sFilename, int& errCode)
{
	std::cout << "\nsize: " << sd->getSize() << "\n";
	FILE *pf_file;
	pf_file = fopen(sFilename.c_str(), "r");
	char name[256];
	int iNum;
	double dNum = 0;
	int blad;
	if (pf_file) {
		//D
		fscanf(pf_file, "%s", name);
		fscanf(pf_file, "%i", &iNum);
		vSetD(iNum, blad);
		std::cout << "D: " << D << "\n";


		//F
		fscanf(pf_file, "%s", name);
		fscanf(pf_file, "%i", &iNum);
		vSetF(iNum, blad);
		std::cout << "F: " << F << "\n";


		//M
		fscanf(pf_file, "%s", name);
		fscanf(pf_file, "%i", &iNum);
		vSetM(iNum, blad);
		std::cout << "M: " << M << "\n";

		//S
		fscanf(pf_file, "%s", name);

		fscanf(pf_file, "%i", &iNum);
		vSetS(iNum, blad);
		std::cout << "S: " << S << "\n";

		//sd
		fscanf(pf_file, "%s", name);
		for (int i = 0; i < D; i++) {
			fscanf(pf_file, "%lf", &dNum);
			sd->setAt(i, dNum, errCode);
		}
		for (int i = 0; i < D; i++) {
			std::cout << "sd[" << i << "]: " << sd->getValueAt(i, errCode) << "\n";
		}
		//sf
		fscanf(pf_file, "%s", name);
		for (int i = 0; i < F; i++) {
			fscanf(pf_file, "%lf", &dNum);
			sf->setAt(i, dNum, errCode);
		}
		for (int i = 0; i < F; i++) {
			std::cout << "sf[" << i << "]: " << sf->getValueAt(i, errCode) << "\n";
		}

		//sm
		fscanf(pf_file, "%s", name);
		for (int i = 0; i < M; i++) {
			fscanf(pf_file, "%lf", &dNum);
			sm->setAt(i, dNum, errCode);

		}
		for (int i = 0; i < M; i++) {
			std::cout << "sm[" << i << "]: " << sm->getValueAt(i, errCode) << "\n";
		}

		//ss
		fscanf(pf_file, "%s", name);
		for (int i = 0; i < S; i++) {
			fscanf(pf_file, "%lf", &dNum);
			ss->setAt(i, dNum, errCode);
		}

		for (int i = 0; i < S; i++) {
			std::cout << "ss[" << i << "]: " << ss->getValueAt(i, errCode) << "\n";
		}

		//cd
		fscanf(pf_file, "%s", name);
		for (int i = 0; i < D; i++) {
			for (int j = 0; j < F; j++) {
				fscanf(pf_file, "%lf", &dNum);
				cd->setAt(i, j, dNum, errCode);
			}
		}

		for (int i = 0; i < D; i++) {
			for (int j = 0; j < F; j++) {
				std::cout << "cd[" << i << "][" << j << "]: " << cd->getValueAt(i, j, errCode) << "\n";
			}
		}

		//cf
		fscanf(pf_file, "%s", name);
		for (int i = 0; i < F; i++) {
			for (int j = 0; j < M; j++) {
				fscanf(pf_file, "%lf", &dNum);
				cf->setAt(i, j, dNum, errCode);
			}
		}

		for (int i = 0; i < F; i++) {
			for (int j = 0; j < M; j++) {
				std::cout << "cf[" << i << "][" << j << "]: " << cf->getValueAt(i, j, errCode) << "\n";
			}
		}

		//cm
		fscanf(pf_file, "%s", name);
		for (int i = 0; i < M; i++) {
			for (int j = 0; j < S; j++) {
				fscanf(pf_file, "%lf", &dNum);
				cm->setAt(i, j, dNum, errCode);
			}
		}

		for (int i = 0; i < M; i++) {
			for (int j = 0; j < S; j++) {
				std::cout << "cm[" << i << "][" << j << "]: " << cm->getValueAt(i, j, errCode) << "\n";
			}
		}


		//ud
		fscanf(pf_file, "%s", name);
		for (int i = 0; i < D; i++) {
			fscanf(pf_file, "%lf", &dNum);
			ud->setAt(i, dNum, errCode);
		}
		for (int i = 0; i < D; i++) {
			std::cout << "ud[" << i << "]: " << ud->getValueAt(i, errCode) << "\n";
		}

		//uf
		fscanf(pf_file, "%s", name);
		for (int i = 0; i < F; i++) {
			fscanf(pf_file, "%lf", &dNum);
			uf->setAt(i, dNum, errCode);
		}
		for (int i = 0; i < F; i++) {
			std::cout << "uf[" << i << "]: " << uf->getValueAt(i, errCode) << "\n";
		}

		//um
		fscanf(pf_file, "%s", name);
		for (int i = 0; i < M; i++) {
			fscanf(pf_file, "%lf", &dNum);
			um->setAt(i, dNum, errCode);
		}
		for (int i = 0; i < M; i++) {
			std::cout << "um[" << i << "]: " << um->getValueAt(i, errCode) << "\n";
		}

		//ps
		fscanf(pf_file, "%s", name);
		for (int i = 0; i < S; i++) {
			fscanf(pf_file, "%lf", &dNum);
			ps->setAt(i, dNum, errCode);
		}
		for (int i = 0; i < S; i++) {
			std::cout << "ps[" << i << "]: " << ps->getValueAt(i, errCode) << "\n";
		}


		//xdminmax

		fscanf(pf_file, "%ls", name);

		for (int i = 0; i < D; i++) {
			for (int j = 0; j < 2 * F; j++) {
				fscanf(pf_file, "%lf", &dNum);
				xdminmax->setAt(i, j, dNum, errCode);

			}
		}

		for (int i = 0; i < D; i++) {
			for (int j = 0; j < 2 * F; j++) {
				std::cout << "xdminmax[" << i << "][" << j << "]: " << xdminmax->getValueAt(i, j, errCode) << "\n";
			}
		}

		//xfminmax

		fscanf(pf_file, "%ls", name);

		for (int i = 0; i < F; i++) {
			for (int j = 0; j < 2 * M; j++) {
				fscanf(pf_file, "%lf", &dNum);
				xfminmax->setAt(i, j, dNum, errCode);

			}
		}

		for (int i = 0; i < F; i++) {
			for (int j = 0; j < 2 * M; j++) {
				std::cout << "xfminmax[" << i << "][" << j << "]: " << xfminmax->getValueAt(i, j, errCode) << "\n";
			}
		}

		//xfminmax

		fscanf(pf_file, "%ls", name);

		for (int i = 0; i < M; i++) {
			for (int j = 0; j < 2 * S; j++) {
				fscanf(pf_file, "%lf", &dNum);
				xmminmax->setAt(i, j, dNum, errCode);
			}
		}

		for (int i = 0; i < M; i++) {
			for (int j = 0; j < 2 * S; j++) {
				std::cout << "xmminmax[" << i << "][" << j << "]: " << xmminmax->getValueAt(i, j, errCode) << "\n";
			}
		}
		fclose(pf_file);

	}

}

void CMscnProblem::saveSolutionFile(std::string sFilename, int& errCode)
{
	FILE* pf_file = fopen(sFilename.c_str(), "w");
	if (pf_file == NULL) {
		errCode = FILE_IS_NULL;
	}

	if (pf_file != NULL) {
		int index = 0;
		//D F M S
		fprintf(pf_file, "D %i", D);
		fprintf(pf_file, "\nF %i", F);
		fprintf(pf_file, "\nM %i", M);
		fprintf(pf_file, "\nS %i", S);

		//xd
		fprintf(pf_file, "\nxd");
		for (int i = 0; i < D; i++) {
			fprintf(pf_file, "\n");
			for (int j = 0; j < F; j++) {
				//fprintf(pf_file, "%f ", xd->getValueAt(i, j, errCode));
				fprintf(pf_file, "%f ", solution->getValueAt(index, errCode));
				index++;
			}
		}

		//xf
		fprintf(pf_file, "\nxf");
		for (int i = 0; i < F; i++) {
			fprintf(pf_file, "\n");
			for (int j = 0; j < M; j++) {
				//fprintf(pf_file, "%f ", xf->getValueAt(i, j, errCode));
				fprintf(pf_file, "%f ", solution->getValueAt(index, errCode));
				index++;

			}
		}

		//xm
		fprintf(pf_file, "\nxm");
		for (int i = 0; i < M; i++) {
			fprintf(pf_file, "\n");
			for (int j = 0; j < S; j++) {
				//fprintf(pf_file, "%f ", xm->getValueAt(i, j, errCode));
				fprintf(pf_file, "%f ", solution->getValueAt(index, errCode));
				index++;

			}
		}

		fclose(pf_file);
	}
}

void CMscnProblem::saveProblemFile(std::string sFilename, int& errCode)
{
	FILE* pf_file = fopen(sFilename.c_str(), "w");
	if (pf_file == NULL) {
		errCode = FILE_IS_NULL;
	}

	if (pf_file != NULL) {

		//D F M S
		fprintf(pf_file, "D %i", D);
		fprintf(pf_file, "\nF %i", F);
		fprintf(pf_file, "\nM %i", M);
		fprintf(pf_file, "\nS %i", S);

		//sd
		fprintf(pf_file, "\nsd\n");
		for (int i = 0; i < D; i++) {
			fprintf(pf_file, "%f ", sd->getValueAt(i, errCode));
		}

		//sf
		fprintf(pf_file, "\nsf\n");
		for (int i = 0; i < F; i++) {
			fprintf(pf_file, "%f ", sf->getValueAt(i, errCode));
		}

		//sm
		fprintf(pf_file, "\nsm\n");
		for (int i = 0; i < M; i++) {
			fprintf(pf_file, "%f ", sm->getValueAt(i, errCode));
		}

		//ss
		fprintf(pf_file, "\nss\n");
		for (int i = 0; i < S; i++) {
			fprintf(pf_file, "%f ", ss->getValueAt(i, errCode));
		}

		//cd
		fprintf(pf_file, "\ncd");
		for (int i = 0; i < D; i++) {
			fprintf(pf_file, "\n");
			for (int j = 0; j < F; j++) {
				fprintf(pf_file, "%f ", cd->getValueAt(i, j, errCode));
			}
		}

		//cf
		fprintf(pf_file, "\ncf");
		for (int i = 0; i < F; i++) {
			fprintf(pf_file, "\n");
			for (int j = 0; j < M; j++) {
				fprintf(pf_file, "%f ", cf->getValueAt(i, j, errCode));
			}
		}

		//cm
		fprintf(pf_file, "\ncm");
		for (int i = 0; i < M; i++) {
			fprintf(pf_file, "\n");
			for (int j = 0; j < S; j++) {
				fprintf(pf_file, "%f ", cm->getValueAt(i, j, errCode));
			}
		}

		//ud
		fprintf(pf_file, "\nud\n");
		for (int i = 0; i < D; i++) {
			fprintf(pf_file, "%f ", ud->getValueAt(i, errCode));
		}
		//uf
		fprintf(pf_file, "\nuf\n");
		for (int i = 0; i < F; i++) {
			fprintf(pf_file, "%f ", uf->getValueAt(i, errCode));
		}
		//um
		fprintf(pf_file, "\num\n");
		for (int i = 0; i < M; i++) {
			fprintf(pf_file, "%f ", um->getValueAt(i, errCode));
		}

		//ps
		fprintf(pf_file, "\np\n");
		for (int i = 0; i < S; i++) {
			fprintf(pf_file, "%f ", ps->getValueAt(i, errCode));
		}

		//xdminMax
		fprintf(pf_file, "\nxdminmax");
		for (int i = 0; i < D; i++) {
			fprintf(pf_file, "\n");
			for (int j = 0; j < 2 * F; j++) {
				fprintf(pf_file, "%f ", xdminmax->getValueAt(i, j, errCode));
			}
		}

		//xfminMax
		fprintf(pf_file, "\nxfminmax");
		for (int i = 0; i < F; i++) {
			fprintf(pf_file, "\n");
			for (int j = 0; j < 2 * M; j++) {
				fprintf(pf_file, "%f ", xfminmax->getValueAt(i, j, errCode));
			}
		}
		//xmminMax
		fprintf(pf_file, "\nxmminmax");
		for (int i = 0; i < M; i++) {
			fprintf(pf_file, "\n");
			for (int j = 0; j < 2 * S; j++) {
				fprintf(pf_file, "%f ", xmminmax->getValueAt(i, j, errCode));
			}
		}
		fclose(pf_file);
	}
	else {
		errCode = FILE_IS_NULL;
	}

}

void CMscnProblem::readSolutionFile(std::string sFilename, int & errCode)
{
		FILE *pf_file;
		pf_file = fopen(sFilename.c_str(), "r");
		char name[256];
		int iNum;
		double dNum = 0;
		int blad;
		int index = 0;
		if (pf_file) {
			//D
			fscanf(pf_file, "%s", name);
			fscanf(pf_file, "%i", &iNum);
			vSetD(iNum, blad);
			std::cout << "D: " << D << "\n";


			//F
			fscanf(pf_file, "%s", name);
			fscanf(pf_file, "%i", &iNum);
			vSetF(iNum, blad);
			std::cout << "F: " << F << "\n";


			//M
			fscanf(pf_file, "%s", name);
			fscanf(pf_file, "%i", &iNum);
			vSetM(iNum, blad);
			std::cout << "M: " << M << "\n";

			//S
			fscanf(pf_file, "%s", name);

			fscanf(pf_file, "%i", &iNum);
			vSetS(iNum, blad);
			std::cout << "S: " << S << "\n";

			//xd
			fscanf(pf_file, "%s", name);
			for (int i = 0; i < D; i++) {
				for (int j = 0; j < F; j++) {
					fscanf(pf_file, "%lf", &dNum);
					xd->setAt(i, j, dNum, errCode);
					solution->setAt(index, dNum, errCode);
					index++;
				}
			}
			for (int i = 0; i < D; i++) {
				for (int j = 0; j < F; j++) {
					std::cout << "xd[" << i << "][" << j << "]: " << xd->getValueAt(i, j, errCode) << "\n";
				}
			}

			//xf
			fscanf(pf_file, "%s", name);
			for (int i = 0; i < F; i++) {
				for (int j = 0; j < M; j++) {
					fscanf(pf_file, "%lf", &dNum);
					xf->setAt(i, j, dNum, errCode);
					solution->setAt(index, dNum, errCode);
					index++;
				}
			}
			for (int i = 0; i < F; i++) {
				for (int j = 0; j < M; j++) {
					std::cout << "xf[" << i << "][" << j << "]: " << xf->getValueAt(i, j, errCode) << "\n";
				}
			}

			//xm
			fscanf(pf_file, "%s", name);
			for (int i = 0; i < M; i++) {
				for (int j = 0; j < S; j++) {
					fscanf(pf_file, "%lf", &dNum);
					xm->setAt(i, j, dNum, errCode);
					solution->setAt(index, dNum, errCode);
					index++;
				}
			}
			for (int i = 0; i < M; i++) {
				for (int j = 0; j < S; j++) {
					std::cout << "xm[" << i << "][" << j << "]: " << xm->getValueAt(i, j, errCode) << "\n";
				}
			}
			fclose(pf_file);
		}
}

int CMscnProblem::getSolutionSize()
{
	return D * F + F * M + M * S;
}

void CMscnProblem::printProblem()
{
	std::cout << "\nD: \n";
	std::cout << D;
	std::cout << "\nF: \n";
	std::cout << F;
	std::cout << "\nM: \n";
	std::cout << M;
	std::cout << "\nS: \n";
	std::cout << S;
	std::cout << "\nsd: \n";
	sd->printArray();
	std::cout << "\nsf: \n";
	sf->printArray();
	std::cout << "\nsm: \n";
	sm->printArray();
	std::cout << "\nss: \n";
	ss->printArray();
	std::cout << "\ncd: \n";
	cd->printMatrix();
	std::cout << "\ncf: \n";
	cf->printMatrix();
	std::cout << "\ncm: \n";
	cm->printMatrix();
	std::cout << "\nud: \n";
	ud->printArray();
	std::cout << "\nuf: \n";
	uf->printArray();
	std::cout << "\num: \n";
	um->printArray();
	std::cout << "\np: \n";
	ps->printArray();
	std::cout << "\nxdminmax: \n";
	xdminmax->printMatrix();
	std::cout << "\nxfminmax: \n";
	xfminmax->printMatrix();
	std::cout << "\nxmminmax: \n";
	xmminmax->printMatrix();
}

void CMscnProblem::vGenerateInstance(int iInstanceSeed, int & errCode)
{
	CRandom crandom(iInstanceSeed);


	for (int i = 0; i < cd->getSizeX(); i++) {
		for (int j = 0; j < cd->getSizeY(); j++) {
			cd->setAt(i, j, crandom.getRandomDouble(MIN_COST, MAX_COST), errCode);
		}
	}

	for (int i = 0; i < cf->getSizeX(); i++) {
		for (int j = 0; j < cf->getSizeY(); j++) {
			cf->setAt(i, j, crandom.getRandomDouble(MIN_COST, MAX_COST), errCode);
		}
	}
	for (int i = 0; i < cm->getSizeX(); i++) {
		for (int j = 0; j < cm->getSizeY(); j++) {
			cm->setAt(i, j, crandom.getRandomDouble(MIN_COST, MAX_COST), errCode);
		}
	}

	for (int i = 0; i < ud->getSize(); i++) {
		ud->setAt(i, crandom.getRandomDouble(MIN_U, MAX_U), errCode);
	}

	for (int i = 0; i < uf->getSize(); i++) {
		uf->setAt(i, crandom.getRandomDouble(MIN_U, MAX_U), errCode);
	}
	for (int i = 0; i < um->getSize(); i++) {
		um->setAt(i, crandom.getRandomDouble(MIN_U, MAX_U), errCode);
	}


	for (int i = 0; i < sd->getSize(); i++) {
		sd->setAt(i, crandom.getRandomDouble(MIN_S, MAX_S), errCode);
	}

	for (int i = 0; i < sf->getSize(); i++) {
		sf->setAt(i, crandom.getRandomDouble(MIN_S, MAX_S), errCode);
	}

	for (int i = 0; i < sm->getSize(); i++) {
		sm->setAt(i, crandom.getRandomDouble(MIN_S, MAX_S), errCode);

	}

	for (int i = 0; i < ss->getSize(); i++) {
		ss->setAt(i, crandom.getRandomDouble(MIN_S, MAX_S), errCode);

	}

	for (int i = 0; i < ps->getSize(); i++) {
		ps->setAt(i, crandom.getRandomDouble(MIN_P, MAX_P), errCode);
	}

	for (int i = 0; i < xdminmax->getSizeX(); i++) {
		for (int j = 0; j < xdminmax->getSizeY(); j += 2) {
			xdminmax->setAt(i, j, crandom.getRandomDouble(MIN_MIN, MAX_MIN), errCode);
			xdminmax->setAt(i, j + 1, crandom.getRandomDouble(MIN_MAX, MAX_MAX), errCode);
		}
	}

	for (int i = 0; i < xfminmax->getSizeX(); i++) {
		for (int j = 0; j < xfminmax->getSizeY(); j += 2) {
			xfminmax->setAt(i, j, crandom.getRandomDouble(MIN_MIN, MAX_MIN), errCode);
			xfminmax->setAt(i, j + 1, crandom.getRandomDouble(MIN_MAX, MAX_MAX), errCode);
		}
	}

	for (int i = 0; i < xmminmax->getSizeX(); i++) {
		for (int j = 0; j < xmminmax->getSizeY(); j += 2) {
			xmminmax->setAt(i, j, crandom.getRandomDouble(MIN_MIN, MAX_MIN), errCode);

			xmminmax->setAt(i, j + 1, crandom.getRandomDouble(MIN_MAX, MAX_MAX), errCode);
		}
	}

}
Array<double>* CMscnProblem::generateSolution(int & errCode)
{
	CRandom crandom;

	int solutionSize = D * F + F * M + M * S;
	Array<double> *solution;
	solution = new Array<double>(solutionSize);

	double sumD = 0;
	double sumF = 0;
	int position = 0;
	for (int i = 0; i < D; i++) {
		for (int j = 0; j < F; j++) {
			int min = 0;
			int max = sd->getValueAt(i, errCode) / (F*0.6);
			double val = crandom.getRandomDouble(min, max);
			solution->setAt(position, val, errCode);
			xd->setAt(i, j, val, errCode);

			position++;
		}
	}


	for (int i = 0; i < F; i++) {
		for (int j = 0; j < M; j++) {
			double xdSumD = 0;

			for (int k = 0; k < D; k++) {
				xdSumD += xd->getValueAt(k, i, errCode);
			}
			xdSumD = xdSumD / (M*0.6);

			int min = 0;
			int max = sf->getValueAt(i, errCode) / (M*0.6);
			if (xdSumD < max) {
				max = xdSumD;
			}
			double val = crandom.getRandomDouble(min, max);
			solution->setAt(position, val, errCode);
			xf->setAt(i, j, val, errCode);

			position++;
			xdSumD = 0;
		}
	}

	for (int i = 0; i < M; i++) {
		for (int j = 0; j < S; j++) {

			double xfSumF = 0;

			for (int k = 0; k < F; k++) {
				xfSumF += xf->getValueAt(k, i, errCode);

			}
			xfSumF = xfSumF /( S*0.6);
			int min = 0;
			int max = sm->getValueAt(i, errCode) / (S*0.6);

			if (ss->getValueAt(j, errCode) / M < max) {
				max = ss->getValueAt(j, errCode) / (M*0.6);

			}
			if (xfSumF < max) {
				max = xfSumF;
			}
			double val = crandom.getRandomDouble(min, max);
			solution->setAt(position, val, errCode);
			xm->setAt(i, j, val, errCode);
			position++;
			xfSumF = 0;
		}
	}

	return solution;
}



Array<double>* CMscnProblem::getSolution(int & errCode)
{
	return solution;
}
void CMscnProblem::setSolutionArray(Array<double>* newSolution,int & errCode)
{
	if (newSolution->getSize() != D * F + F * M + M * S) {
		errCode = WRONG_SIZE;
	}
	else {
		solution = newSolution;
	}
}
void CMscnProblem::setSolution(double* newSolution)
{
	int errCode = 0;
	for (int i = 0; i < D*F + F * M + M * S; i++) {
		solution->setAt(i, newSolution[i], errCode);

	}
}

/*
bool CMscnProblem::readFromFileSemicolon(std::string sFilename)
{
	FILE *pf_file;
	pf_file = fopen(sFilename.c_str(), "r");
	char name[256];
	char semicolon;
	int iNum = 0;
	double dNum = 0;

	int errCode = 0;
	if (pf_file) {
		//D
		fscanf(pf_file, "%s", name);
		fscanf(pf_file, "%i", &iNum);
		vSetD(iNum, errCode);
		std::cout << "D: " << D << "\n";


		//F
		fscanf(pf_file, "%s", name);
		fscanf(pf_file, "%i", &iNum);
		vSetF(iNum, errCode);
		std::cout << "F: " << F << "\n";


		//M
		fscanf(pf_file, "%s", name);
		fscanf(pf_file, "%i", &iNum);
		vSetM(iNum, errCode);
		std::cout << "M: " << M << "\n";

		//S
		fscanf(pf_file, "%s", name);
		fscanf(pf_file, "%i", &iNum);
		vSetS(iNum, errCode);
		std::cout << "S: " << S << "\n";


		//sd
		fscanf(pf_file, "%s", name);

		for (int i = 0; i < D; i++) {
			fscanf(pf_file, "%lf", &dNum);
			vSetSd(dNum, i);
			fscanf(pf_file, "%c", &semicolon);
			if (semicolon != ';') {
				return false;
			}

		}
		for (int i = 0; i < D; i++) {
			std::cout << "sd[" << i << "]: " << sd[i] << "\n";
		}


		//sf
		fscanf(pf_file, "%s", name);
		for (int i = 0; i < F; i++) {
			fscanf(pf_file, "%lf", &dNum);
			vSetSf(dNum, i);
			fscanf(pf_file, "%c", &semicolon);
			if (semicolon != ';') {
				return false;
			}
		}
		for (int i = 0; i < F; i++) {
			std::cout << "sf[" << i << "]: " << sf[i] << "\n";
		}

		//sm
		fscanf(pf_file, "%s", name);
		for (int i = 0; i < M; i++) {
			fscanf(pf_file, "%lf", &dNum);
			vSetSm(dNum, i);
			fscanf(pf_file, "%c", &semicolon);
			if (semicolon != ';') {
				return false;
			}
		}
		for (int i = 0; i < M; i++) {
			std::cout << "sm[" << i << "]: " << sm[i] << "\n";
		}

		//ss
		fscanf(pf_file, "%s", name);
		for (int i = 0; i < S; i++) {
			fscanf(pf_file, "%lf", &dNum);
			vSetSs(dNum, i);
			fscanf(pf_file, "%c", &semicolon);
			if (semicolon != ';') {
				return false;
			}
		}

		for (int i = 0; i < S; i++) {
			std::cout << "ss[" << i << "]: " << ss[i] << "\n";
		}

		//cd
		fscanf(pf_file, "%s", name);
		for (int i = 0; i < D; i++) {
			for (int j = 0; j < F; j++) {
				fscanf(pf_file, "%lf", &dNum);
				vSetCd(dNum, i, j);
				fscanf(pf_file, "%c", &semicolon);
				if (semicolon != ';') {
					return false;
				}
			}
		}

		for (int i = 0; i < D; i++) {
			for (int j = 0; j < F; j++) {
				std::cout << "cd[" << i << "][" << j << "]: " << cd[i][j] << "\n";
			}
		}

		//cf
		fscanf(pf_file, "%s", name);
		for (int i = 0; i < F; i++) {
			for (int j = 0; j < M; j++) {
				fscanf(pf_file, "%lf", &dNum);
				vSetCf(dNum, i, j);
				fscanf(pf_file, "%c", &semicolon);
				if (semicolon != ';') {
					return false;
				}
			}
		}

		for (int i = 0; i < F; i++) {
			for (int j = 0; j < M; j++) {
				std::cout << "cf[" << i << "][" << j << "]: " << cf[i][j] << "\n";
			}
		}

		//cm
		fscanf(pf_file, "%s", name);
		for (int i = 0; i < M; i++) {
			for (int j = 0; j < S; j++) {
				fscanf(pf_file, "%lf", &dNum);
				vSetCm(dNum, i, j);
				fscanf(pf_file, "%c", &semicolon);
				if (semicolon != ';') {
					return false;
				}
			}
		}

		for (int i = 0; i < M; i++) {
			for (int j = 0; j < S; j++) {
				std::cout << "cm[" << i << "][" << j << "]: " << cm[i][j] << "\n";
			}
		}


		//ud
		fscanf(pf_file, "%s", name);
		for (int i = 0; i < D; i++) {
			fscanf(pf_file, "%lf", &dNum);
			//vSetUd(dNum, i);
			ud[i] = dNum;
			fscanf(pf_file, "%c", &semicolon);
			if (semicolon != ';') {
				return false;
			}
		}
		for (int i = 0; i < D; i++) {
			std::cout << "ud[" << i << "]: " << ud[i] << "\n";
		}

		//uf
		fscanf(pf_file, "%s", name);
		for (int i = 0; i < F; i++) {
			fscanf(pf_file, "%lf", &dNum);
			//vSetUf(dNum, i);
			uf[i] = dNum;
			fscanf(pf_file, "%c", &semicolon);
			if (semicolon != ';') {
				return false;
			}
		}
		for (int i = 0; i < F; i++) {
			std::cout << "uf[" << i << "]: " << uf[i] << "\n";
		}

		//um
		fscanf(pf_file, "%s", name);
		for (int i = 0; i < M; i++) {
			fscanf(pf_file, "%lf", &dNum);
			//vSetUm(dNum, i);
			um[i] = dNum;
			fscanf(pf_file, "%c", &semicolon);
			if (semicolon != ';') {
				return false;
			}
		}
		for (int i = 0; i < M; i++) {
			std::cout << "um[" << i << "]: " << um[i] << "\n";
		}

		//ps
		fscanf(pf_file, "%s", name);
		for (int i = 0; i < S; i++) {
			fscanf(pf_file, "%lf", &dNum);
			//vSetPs(dNum, i);
			ps[i] = dNum;
			fscanf(pf_file, "%c", &semicolon);
			if (semicolon != ';') {
				return false;
			}
		}
		for (int i = 0; i < S; i++) {
			std::cout << "ps[" << i << "]: " << ps[i] << "\n";
		}


		//xdminmax

		fscanf(pf_file, "%ls", name);

		for (int i = 0; i < D; i++) {
			for (int j = 0; j < 2 * F; j++) {
				fscanf(pf_file, "%lf", &dNum);
				xdminmax[i][j] = dNum;
				fscanf(pf_file, "%c", &semicolon);
				if (semicolon != ';') {
					return false;
				}
			}

		}

		for (int i = 0; i < D; i++) {
			for (int j = 0; j < 2 * F; j++) {
				std::cout << "xdminmax[" << i << "][" << j << "]: " << xdminmax[i][j] << "\n";
			}
		}

		//xfminmax

		fscanf(pf_file, "%ls", name);

		for (int i = 0; i < F; i++) {
			for (int j = 0; j < 2 * M; j++) {
				fscanf(pf_file, "%lf", &dNum);
				xfminmax[i][j] = dNum;
				fscanf(pf_file, "%c", &semicolon);
				if (semicolon != ';') {
					return false;
				}
			}

		}

		for (int i = 0; i < F; i++) {
			for (int j = 0; j < 2 * M; j++) {
				std::cout << "xfminmax[" << i << "][" << j << "]: " << xfminmax[i][j] << "\n";
			}
		}

		//xfminmax

		fscanf(pf_file, "%ls", name);

		for (int i = 0; i < M; i++) {
			for (int j = 0; j < 2 * S; j++) {
				fscanf(pf_file, "%lf", &dNum);
				xmminmax[i][j] = dNum;
				fscanf(pf_file, "%c", &semicolon);
				if (semicolon != ';') {
					return false;
				}
			}

		}

		for (int i = 0; i < M; i++) {
			for (int j = 0; j < 2 * S; j++) {
				std::cout << "xmminmax[" << i << "][" << j << "]: " << xmminmax[i][j] << "\n";
			}
		}

		fclose(pf_file);
	}
}
*/