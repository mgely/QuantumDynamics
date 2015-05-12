#include "header.h"
using namespace std;

int main()
{	

	int N_x = 256;        			// number of x-axis grid points
	int N_y = N_x;      
	double L_x = 100;              	// system exts from x=0 to x=L
	double L_y = 100;              
	int initial_conditions = TUNN; 	// Choose from FREE, HARM (harmonic), TUNN (tunneling)
	int output = NORM;				// Choose from NORM, REAL, IMAG for plotting 
	double dt = 0.1;
	double simulation_time = 100;

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