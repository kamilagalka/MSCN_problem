#include "pch.h"
#include "CRandom.h"
#include <iostream>
#include <random>


CRandom::CRandom()
{
	gen = new std::mt19937(rd());
}

CRandom::CRandom(int seed)
{
	gen = new std::mt19937(seed);
}

CRandom::~CRandom()
{
	delete gen;
}

int CRandom::getRandomInt(int lower_bound, int upper_bound)
{
	std::uniform_int_distribution<> dis(lower_bound, upper_bound);

	return dis(*gen);
}

double CRandom::getRandomDouble(double lower_bound, double upper_bound)
{
	std::uniform_real_distribution<> dis(lower_bound, upper_bound);
	return dis(*gen);
}