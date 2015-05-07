#include "header.h"

void quantum::step(){
	FFT FFT_instance = FFT(N);
	cpxnumber temp;
	double k;
	wave = FFT_instance.transform(wave);

	for (int i = 1; i < N/2; ++i){
		k = 0.5*(double)N/(double)i;
		temp = cpxnumber(-k*k*dt/4.);// Not right ?
		wave[i].mult_cpx(temp);
		wave[N-i-1].mult_cpx(temp);
	}
	
	wave = FFT_instance.itransform(wave);
}

void quantum::simulation(std::ostream& file,int time_steps){
	for (int i = 0; i < time_steps; ++i){
		for (int j = 0; j < N; ++j)
			file<<wave[j].norm()<<" ";
		file<<"\n";
		step();
	}
}