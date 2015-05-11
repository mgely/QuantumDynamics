#include "header.h"
using namespace std;

simulation::simulation(int N_x, int N_y, double L_x, double L_y, double dt, double simulation_time, int initial_conditions, int output) :
N_x(N_x), N_y(N_y),L_x(L_x), L_y(L_y), dt(dt), simulation_time(simulation_time), initial_conditions(initial_conditions), output(output)
{	
	cout << "INITIALIZING..." << endl;
	set_patial_grid();
	allocate_memory();

	setV();
	setT();
	setpsi();


	cout << "OPTIMIZING FFTW..." << endl;
	// FFT choice of plan
	fft = fftw_plan_dft_2d(N_x, N_y,psi,psi,FFTW_FORWARD, FFTW_PATIENT);
	setpsi();

	ifft = fftw_plan_dft_2d(N_x, N_y,psi,psi,FFTW_BACKWARD, FFTW_PATIENT);
	setpsi();
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
			psi[i*N_x + j][1] = (psi[i*N_x + j][0] * V_exp_factor[i*N_x + j][1] + psi[i*N_x + j][1] * V_exp_factor[i*N_x + j][0])/((double)N_x*(double)N_y);
		}
	}

//        	for (int i = 0; i < N_x; ++i)
// {
// 	for (int j = 0; j < N_y; ++j)
// 	{
// 		cout << psi[i*N_x+j][0] << '\t';
// 	}
// 	cout << endl;
// }

}

void simulation::print_to_file(int flag){
	FILE* to_plot = fopen("plot.dat","w");

	if(flag == 0){ // "norm" case (returns psi squared)
		for (int i = 0; i < N_x; ++i){
			for (int j = 0; j < N_y; ++j)
				fprintf(to_plot, "%f %f %f\n", x[i] , y[j] , psi[i*N_x + j][0]*psi[i*N_x + j][0] + psi[i*N_x + j][1]*psi[i*N_x + j][1]);
			fprintf(to_plot,"\n");
		}
	}

	if(flag == 1){ // "real" case
		for (int i = 0; i < N_x; ++i){
			for (int j = 0; j < N_y; ++j)
				fprintf(to_plot, "%f %f %f\n", x[i] , y[j] , psi[i*N_x + j][0]);
			fprintf(to_plot,"\n");
		}
	}

	if(flag == 2){ // "imag" case
		for (int i = 0; i < N_x; ++i){
			for (int j = 0; j < N_y; ++j)
				fprintf(to_plot, "%f %f %f\n", x[i] , y[j] , psi[i*N_x + j][1]);
			fprintf(to_plot,"\n");
		}
	}

	fclose(to_plot);
}

void simulation::simulate(){
	
	double t = 0;
	pid_t gnupid = fork(); // After fork() has been called there are two different processes running:
							// the parent process identified by gnupid > 0, and the child process gnupid == 0

	// Child process : gnuplot
	if (gnupid == 0)
	{
		ofstream cmdfile("gnusettings2d.txt");

		cmdfile   << "set term x11" << "\n"
				  << "set title \" time dependent Schrodinger wavefunction\" " << "\n"
				  << "set time" << "\n"
				  << "set xtics" <<"\n"
				  << "set ytics" << "\n"
				  << "set xrange [0:" << L_x << "]" << "\n"
				  << "set yrange [0:" << L_y << "]" << "\n"
				  << "plot \"plot.dat\" using 1:2:3 with image" << "\n"
				  << "pause 0.001 " << "\n"
				  << "reread;" ;

		    cmdfile.flush();
			cmdfile.close();

			execlp("gnuplot","gnuplot","gnusettings2d.txt", (char*) NULL);
	}

	// Parent process: the solver
	else if (gnupid > 0)
	{
		mkfifo("plot.dat", S_IWUSR | S_IRUSR); // making a FIFO file to write data for gnuplot
		while (t<simulation_time)
		{	
			cout << "time : "<< t << endl;
			step();
			print_to_file(NORM);
			t += dt;
		}
	}

	else if ( gnupid == -1)
		fprintf(stdout," Error forking process");

	kill(gnupid,SIGTERM);
	cout << "END OF SIMULATION" << endl;
}