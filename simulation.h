#ifndef SIMULATION_H
#define SIMULATION_H

#include <fftw3.h>
const int NORM = 0;
const int REAL = 1;
const int IMAG = 2;

class simulation{
public: 
	int N_x,N_y; 
	fftw_complex * psi;
	fftw_complex * V_exp_factor;
	fftw_complex * T_exp_factor;
	fftw_plan fft;
	fftw_plan ifft;

public:
	simulation();
	~simulation () {};

	void step();
	void print_to_file(int flag);
};

#endif