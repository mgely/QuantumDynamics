#ifndef SIMULATION_H
#define SIMULATION_H

#include <fftw3.h>
#include <vector>



class simulation{
private: 
	int N_x,N_y; 
	fftw_complex * psi;
	fftw_complex * V_exp_factor;
	fftw_complex * T_exp_factor;
	fftw_plan fft;
	fftw_plan ifft;
	std::vector<double> x;
	std::vector<double> y;
	double L_x;             // system exts from x=0 to x=L
	double L_y;             // system exts from y=0 to y=L
	double dt;
	double simulation_time;
	int initial_conditions;
	int output;

public:
	simulation(int N_x, int N_y, double L_x, double L_y, double dt, double simulation_time, int initial_conditions, int output);
	~simulation () {};

	void step();
	void print_to_file(int flag);
	void simulate();
	void setV();
	void setT();
	void setpsi();
	void set_patial_grid();
	void allocate_memory();
};

#endif