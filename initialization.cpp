#include "header.h"
using namespace std;

void simulation::allocate_memory(){
	V_exp_factor = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N_x*N_y);
	T_exp_factor = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N_x*N_y);
	psi = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N_x*N_y);
}

void simulation::setV(){

	double V[N_x][N_y];

	if (initial_conditions == FREE){
		for (int i = 0; i < N_x; ++i){
			for (int j = 0; j < N_y; ++j){
				V[i][j] = 0;
			}
		}
	}

	else if (initial_conditions == HARM){
		double V_x_0 = L_x/2;
    	double V_y_0 = L_y/2;
		for (int i = 0; i < N_x; ++i){
			for (int j = 0; j < N_y; ++j){
				V[i][j] = ((x[i]-V_x_0)*(x[i]-V_x_0)+(y[j]-V_y_0)*(y[j]-V_y_0));
			}
		}
	}

	else if (initial_conditions == TUNN){
		double V_width = 2;
		double V_center = L_x/2+10.;
		double half_width = abs(0.5 * V_width);
		double V_height = 2;
		for (int i = 0; i < N_x; ++i){
			for (int j = 0; j < N_y; ++j)
				V[i][j] = V_height * exp(-(x[i]-V_center)*(x[i]-V_center)/ V_width);
		}
	}  

	else if (initial_conditions == SLIT){
		double V_width = 1;
		double slit_width = 1;
		double V_center = L_x/2.;
		double half_width = abs(0.5 * V_width);
		double V_height = 1000;
		for (int i = 0; i < N_x; ++i){
			for (int j = 0; j < N_y; ++j)
				if(abs(y[i] - L_y/2.)< slit_width)
					V[i][j] = 0;
				else
					V[i][j] = V_height * exp(-(x[i]-V_center)*(x[i]-V_center)/ V_width);
		}
	}  

	for (int i = 0; i < N_x; ++i)
	{
		for (int j = 0; j < N_y; ++j){
			V_exp_factor[i*N_x + j][0] = cos(-V[i][j]/2.*dt);
			V_exp_factor[i*N_x + j][1] = sin(-V[i][j]/2.*dt);
		}
	}
}

void simulation::setT(){

	double px;
	double py;
	double theta;
	for (int i = 0; i < N_x; ++i){
		for (int j = 0; j < N_y; ++j){
			if(i<N_x/2)
				px = i;
			else 
				px = i-N_x;
			px = px*2*PI/L_x;

			if(j<N_y/2)
				py = j;
			else
				py = j-N_y;
			py = py*2*PI/L_y;

			theta = - (px * px + py * py) / 2 * dt;
			T_exp_factor[i*N_x + j][0] = cos(theta);
			T_exp_factor[i*N_x + j][1] = sin(theta);
		}
	}
}

void simulation::setpsi(){

double x_0,y_0,sigma_0,k_0;
double gaussian[N_x][N_y];

	if (initial_conditions == FREE){
		x_0 = L_x/2.; // location of center
		y_0 = L_y/2.;
		sigma_0 = 5.*L_x / 100.; // width of wave packet
		k_0 = 0; // average wavenumber in the x - direction
	}

	if (initial_conditions == HARM){
		x_0 = L_x/2.-10.; // location of center
		y_0 = L_y/2.;
		sigma_0 = 5.*L_x / 100.; // width of wave packet
		k_0 = 0; // average wavenumber in the x - direction
	}

	if (initial_conditions == TUNN || initial_conditions == SLIT){

		x_0 = L_x/2.-10.; // location of center
		y_0 = L_y/2.;
		sigma_0 = 10.*L_x / 100.; // width of wave packet
		k_0 = 1.5; // average wavenumber in the x - direction
	}


	for (int i = 0; i < N_x; ++i)
	{
		for (int j = 0; j < N_y; ++j){
			gaussian[i][j] = exp(-((x[i]-x_0)*(x[i]-x_0)+(y[j]-y_0)*(y[j]-y_0))/ (2 * sigma_0 * sigma_0));
			psi[i*N_x + j][0] = cos(k_0*x[i])*gaussian[i][j];
			psi[i*N_x + j][1] = sin(k_0*x[i])*gaussian[i][j];
		}
	}
}

void simulation::set_patial_grid(){
	x.resize(N_x);
	double h_x = L_x/(double)N_x;
	for (int i = 0; i < N_x; ++i)
	{
		x[i] = ((double)i+1) * h_x;
	}

	y.resize(N_y);
	double h_y = L_y/(double)N_y;
	for (int j = 0; j < N_y; ++j)
	{
		y[j] = ((double)j+1) * h_y;
	}
}