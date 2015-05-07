#include "header.h"
using namespace std;

int main()
{
	int N = 1024; // Size of input function, should always be a 
	ofstream file("data.dat");
	int time_steps = 2;
	double dt = 0.002;

	std::vector<cpxnumber> wave (N);
	for (int i = 0; i < N; ++i)
	 	wave[i] = cpxnumber(exp(-(10.*((double)i-(double)N/2.)/(double)N)*(10.*((double)i-(double)N/2.)/(double)N)),0);
		// wave[i] = cpxnumber(((double)i)/4.);

	quantum q = quantum(wave,N,dt);
	q.simulation(file,time_steps);
	
	return 0;
}