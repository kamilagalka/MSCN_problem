#pragma once
#include <random>
class CRandom
{
public:
	CRandom();
	CRandom(int seed);
	~CRandom();
	int getRandomInt(int lower_bound, int upper_bound);
	double getRandomDouble(double lower_bound, double upper_bound);

private:
	std::random_device rd;
	std::mt19937 *gen;
};
