#ifndef QUANTUM_H
#define QUANTUM_H

#include <cmath>
#include <vector>
#include <iostream>
#include "complex.h"

#define PI 3.14159265358979323846

class quantum{
protected: 
	int N; 	// Number of samples in the FFT
	double dt;
	std::vector<cpxnumber> wave;
	

public:
	quantum(std::vector<cpxnumber> wave,int N, double dt) : wave(wave), N(N), dt(dt) {}
	~quantum () {};

	void step();
	void simulation(std::ostream& file,int time_steps);

};

#endif