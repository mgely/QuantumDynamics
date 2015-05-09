#include "header.h"
using namespace std;

int main()
{	

	// The spatial grid
	int N_x = 256;             // number of x-axis grid points
	int N_y = N_x;             // number of x-axis grid points
	double L_x = 100;              // system exts from x=0 to x=L
	double L_y = 100;              // system exts from y=0 to y=L
	int initial_conditions = FREE;
	int output = NORM;	
	double dt = 0.2;
	double simulation_time = 4;

	simulation s(N_x,
		N_y,
		L_x,
		L_y,
		dt, 
		simulation_time,
		initial_conditions,
		output);

	s.simulate();
	return 0;
}