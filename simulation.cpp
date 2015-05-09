#include "header.h"
using namespace std;

simulation::simulation(){

	ifstream V_exp_factor_file("V_exp_factor.dat");
	ifstream T_exp_factor_file("T_exp_factor.dat");
	ifstream psi_file("psi.dat");
	ifstream dimensions("dimensions.dat");

	// Dimensions
	dimensions >> N_x >> N_y;

	// Memory allocation
	V_exp_factor = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N_x*N_y);
	T_exp_factor = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N_x*N_y);
	psi = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N_x*N_y);



	fftw_complex * psi_temp;
	psi_temp = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N_x*N_y);

	// Psi, V, T setup
	for (int i = 0; i < N_x; ++i){	
		for (int j = 0; j < N_y; ++j){
			V_exp_factor_file >> V_exp_factor[i*N_x + j][0] >> V_exp_factor[i*N_x + j][1];
			T_exp_factor_file >> T_exp_factor[i*N_x + j][0] >> T_exp_factor[i*N_x + j][1];
			psi_file >> psi_temp[i*N_x + j][0] >> psi_temp[i*N_x + j][1];
		}
	}

	// FFT choice of plan
	
	for (int i = 0; i < N_x*N_y; ++i){	
		for (int j = 0; j < 2; ++j)
			psi[i][j] = psi_temp[i][j];
	}
	fft = fftw_plan_dft_2d(N_x, N_y,psi,psi,1, FFTW_PATIENT);

	for (int i = 0; i < N_x*N_y; ++i){	
		for (int j = 0; j < 2; ++j)
			psi[i][j] = psi_temp[i][j];
	}
	ifft = fftw_plan_dft_2d(N_x, N_y,psi,psi,-1, FFTW_PATIENT);

	for (int i = 0; i < N_x*N_y; ++i){	
		for (int j = 0; j < 2; ++j)
			psi[i][j] = psi_temp[i][j];
	}

	// Close files
	V_exp_factor_file.close();
	T_exp_factor_file.close();
	psi_file.close();
	dimensions.close();
}

void simulation::step(){

// first half potential phase rotation

	for (int i = 0; i < N_x; ++i)
	{	
		for (int j = 0; j < N_y; ++j)
		{
			psi[i*N_x + j][0] = psi[i*N_x + j][0] * V_exp_factor[i*N_x + j][0] - psi[i*N_x + j][1] * V_exp_factor[i*N_x + j][1];
			psi[i*N_x + j][1] = psi[i*N_x + j][0] * V_exp_factor[i*N_x + j][1] + psi[i*N_x + j][1] * V_exp_factor[i*N_x + j][0];
		}
	}
    
    
// FFT to momentum space
    fftw_execute(fft);
    
// kinetic phase rotation
	for (int i = 0; i < N_x; ++i)
	{	
		for (int j = 0; j < N_y; ++j)
		{
			psi[i*N_x + j][0] = psi[i*N_x + j][0] * T_exp_factor[i*N_x + j][0] - psi[i*N_x + j][1] * T_exp_factor[i*N_x + j][1];
			psi[i*N_x + j][1] = psi[i*N_x + j][0] * T_exp_factor[i*N_x + j][1] + psi[i*N_x + j][1] * T_exp_factor[i*N_x + j][0];
		}
	}
    
// FFT back to position space
    fftw_execute(ifft);
    
// second half potential phase rotation
	for (int i = 0; i < N_x; ++i)
	{	
		for (int j = 0; j < N_y; ++j)
		{
			psi[i*N_x + j][0] = (psi[i*N_x + j][0] * V_exp_factor[i*N_x + j][0] - psi[i*N_x + j][1] * V_exp_factor[i*N_x + j][1])/((double)N_x*(double)N_y);
			psi[i*N_x + j][1] = (psi[i*	N_x + j][0] * V_exp_factor[i*N_x + j][1] + psi[i*N_x + j][1] * V_exp_factor[i*N_x + j][0])/((double)N_x*(double)N_y);
		}
	}
}

void simulation::print_to_file(int flag){
	ofstream psi_file("psi.dat");

	if(flag == 0){ // "norm" case (returns psi squared)
		for (int i = 0; i < N_x; ++i)
		{	
			for (int j = 0; j < N_y; ++j)
			{	
				psi_file << psi[i*N_x + j][0]*psi[i*N_x + j][0] + psi[i*N_x + j][1]*psi[i*N_x + j][1] << ' ';
			}
			psi_file << '\n';
		}
	}

	if(flag == 1){ // "real" case
		for (int i = 0; i < N_x; ++i)
		{	
			for (int j = 0; j < N_y; ++j)
			{	
				psi_file << psi[i*N_x + j][0] << ' ';
			}
			psi_file << '\n';
		}
	}

	if(flag == 2){ // "imag" case
		for (int i = 0; i < N_x; ++i)
		{	
			for (int j = 0; j < N_y; ++j)
			{	
				psi_file << psi[i*N_x + j][1] << ' ';
			}
			psi_file << '\n';
		}
	}

	psi_file.close();
}

	// for (int i = 0; i < N_x; ++i){	
	// 	for (int j = 0; j < N_y; ++j){
	// 		cout << psi[i*N_x + j][0] << ' ';
	// 	}
	// 	cout << '\n';
	// }