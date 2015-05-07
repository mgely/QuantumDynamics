#include "header.h"

FFT::FFT(int N_in) {

	N = N_in;

	// Compute the twiddle factors
	twiddle.resize(N);

	for (int k = 0; k < N/2; ++k)
	{
		twiddle[k] = cpxnumber(-2.*PI*(double)k/(double)N);
		twiddle[k+N/2] = cpxnumber(2.*PI*(double)k/(double)N);
	}

}

std::vector<cpxnumber> FFT::transform(std::vector<cpxnumber> in, int N_out){
		if(N_out == 0)
			N_out = N;

		std::vector<cpxnumber> out (N_out);

		if(N_out == 1){
			out[0] = in[0];
		}

		else{

			std::vector<cpxnumber> in_pair (N_out/2);
			std::vector<cpxnumber> in_impair (N_out/2);

			for (int i = 0; i < N_out/2 ; ++i){

				in_pair[i] = in[2*i];
				in_impair[i] = in[2*i+1];
			}

			// For speed, I could avoid passing vectors as arguments and use some 
			// kind of argument that is enough to know the return in case N == 0
			std::vector<cpxnumber> out_first  = transform(in_pair,N_out/2);
			std::vector<cpxnumber> out_second = transform(in_impair,N_out/2);

			for (int i = 0; i < N_out/2 ; ++i){

				out_second[i].mult_cpx(twiddle[N*i/N_out]);

				out[i].add(out_first[i],out_second[i]);
				out[i+N_out/2].substract(out_first[i],out_second[i]);
			}
		}
		
	return out;
}

std::vector<cpxnumber> FFT::itransform(std::vector<cpxnumber> in, int N_out){
		if(N_out == 0)
			N_out = N;
		
		std::vector<cpxnumber> out (N_out);

		if(N_out == 1){
			out[0] = in[0];
			out[0].mult_real(1./(double)N);
		}

		else{

			std::vector<cpxnumber> in_pair (N_out/2);
			std::vector<cpxnumber> in_impair (N_out/2);

			for (int i = 0; i < N_out/2 ; ++i){

				in_pair[i] = in[2*i];
				in_impair[i] = in[2*i+1];
			}

			// For speed, I could avoid passing vectors as arguments and use some 
			// kind of argument that is enough to know the return in case N == 0
			std::vector<cpxnumber> out_first  = transform(in_pair,N_out/2);
			std::vector<cpxnumber> out_second = transform(in_impair,N_out/2);

			for (int i = 0; i < N_out/2 ; ++i){

				out_second[i].mult_cpx(twiddle[N-N*i/N_out]);

				out[i].add(out_first[i],out_second[i]);
				out[i+N_out/2].substract(out_first[i],out_second[i]);
			}
		}
		
	return out;
}