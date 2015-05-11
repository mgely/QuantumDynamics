#ifndef HEADER_H
#define HEADER_H

#include <iostream>
#include <unistd.h>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include "simulation.h"
#include <fftw3.h>
#include <signal.h>
#include <sys/stat.h> 

// ouptut flags
const int NORM = 0;
const int REAL = 1;
const int IMAG = 2;

// initial conditions flags
const int FREE = 0;
const int HARM = 1;
const int TUNN = 2;
const int SLIT = 3;
const int DOUB = 4;

const double PI = 3.14159265358979323846;

#endif
