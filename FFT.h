#ifndef FFT_H
#define FFT_H

#include <cmath>
#include <vector>
#include "complex.h"

#define PI 3.14159265358979323846

class FFT{
protected: 
	int N; 							// Number of samples in the FFT
	std::vector<cpxnumber> twiddle;
	

public:
	FFT(int N);
	~FFT () {};

	std::vector<cpxnumber> transform(std::vector<cpxnumber> f_in, int N_out = 0);
	std::vector<cpxnumber> itransform(std::vector<cpxnumber> f_in, int N_out = 0);
};

#endif