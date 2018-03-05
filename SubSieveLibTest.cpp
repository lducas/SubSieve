#include "SubSieveLib.cpp"
#include <stdio.h>
#include <string>
#include <iostream>
#include <fstream>
#include <assert.h>

int main( int argc, const char* argv[] )
{
	assert(argc == 2); 
	int n = atoi( argv[1] );

	string filename = "./svpchallenge/" + to_string(n) + ".gso";
	ifstream gso_file;
	gso_file.open(filename);

	double gh;
	double* mu = new double[n*n];

	gso_file >> gh;
	cerr << "gh" << gh << endl;

	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < n; ++j)
		{
			gso_file >> mu[n*i + j];
			//cerr << mu[n*i + j] << " ";
		}
		//cerr << endl;
	}

	initialize(n, 0, mu, gh);
	sieve(0);

}