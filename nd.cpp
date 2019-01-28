#include "stdafx.h"
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <algorithm>
#include <windows.h>
#include <random>
#include <cmath> // atan2

using namespace std;

// Monte Carlo simulation
static double MinFromN(int index, double* D, int n, int precision)
{
	long long N = (long long)pow(10, precision);
	std::default_random_engine generator;

	int success = 0;
	for (int i = 0; i<N; i++) // 10^precision tests generation
	{
		int minj = -1;
		double min = 1e9; // upper margin for all standard deviations
		double f = 0;

		for (int j = 0; j<n; j++) // one test
		{
			std::normal_distribution<double> distribution(0, D[j]);
			f = distribution(generator); // xj value for the current test
			f *= f; // square the values to make the comparison by modulos
			if (f<min) // keep track of the square of the smallest one ( min ) and it's index ( minj )
			{
				minj = j;
				min = f;
			}
		}
		if (index == minj)
		{
			success++; // count the tests when index came out smallest
		}
	}
	return (double)(success) / N; // return the probability

}

int main()
{
	int index = 0;
	//double D[] = { 9,2,14,18,6,13,19,60,12,30 };
	double D[] = { 9,2,14,18,6,13,19,60,12,30,100,7,16,2,5,65,32,12,8,3 };
	int N = sizeof(D) / sizeof(D[0]);
	double pi = 3.1415926535897932384626433;

	double* match = new double[N];
	int count = 0;
	double sum = 0;
	for (int i = 0; i < N; ++i) // Don Reba approach, utilizing Owen's T-function
	{
		match[count] = 2 * atan2((D[i] / D[0]), 1) / pi;
		fprintf_s(stdout, "match %d %f\n", i, match[count]);
		sum += (1 - match[count]) / (match[count]);
		count++;
	}

	double MonteCarlo = 1;
	MonteCarlo = MinFromN(index, D, N, 6);
	fprintf_s(stdout, "\nMonte Carlo   - >%f\n", MonteCarlo);

	// fprintf_s(Out, "\nMe tripple elements - >%f\n", 1- 2*atan2((D[0]/D[1]),1)/pi- 2*atan2((D[0]/D[2]),1)/pi + 2*atan2(D[0]*D[0]/sqrt(D[0]*D[0]*D[1]*D[1]+D[0]*D[0]*D[2]*D[2]+D[1]*D[1]*D[2]*D[2]),1)/pi);
	// fprintf_s(Out, "Me couple elements - >%f\n",1-2*atan2((D[0]/D[1]),1)/pi);
	// fprintf_s(Out, "Don Reba couple elements - >%f\n",  2*atan2((D[1]/D[0]),1)/pi);

	double Pmagic = 0;
	double Sone = D[0];
	double Stwo = D[1];
	double Sthree = 0;
	for (int i = 2; i < N; ++i) // My Induction, utilizing Owen's T-function and erf-functions properties
	{
		Sthree = D[i];
		double tempPtripple = 1 - 2 * atan2((Sone / Stwo), 1) / pi - 2 * atan2((Sone / Sthree), 1) / pi + 2 * atan2(Sone*Sone / sqrt(Sone*Sone*Stwo*Stwo + Sone * Sone*Sthree*Sthree + Stwo * Stwo*Sthree*Sthree), 1) / pi;
		Stwo = Sone * tan(pi*tempPtripple / 2);
	}

	double induction = 2 * atan2((Stwo / Sone), 1) / pi;
	fprintf_s(stdout, "\nMy  Induction - >%f %f", induction, (induction - MonteCarlo) * 100 / MonteCarlo);
	fprintf_s(stdout, "\nDon Reba      - >%f %f", 1 / sum, (1 / sum - MonteCarlo) * 100 / MonteCarlo);
	fprintf_s(stdout, "\nMix           - >%f %f \n", (induction + 1 / sum) / 2, ((induction + 1 / sum) / 2 - MonteCarlo) * 100 / MonteCarlo);


	delete[] match;

}
