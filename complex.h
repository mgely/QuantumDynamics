#ifndef COMPLEX_H
#define COMPLEX_H

#include <cmath>

class cpxnumber{
protected: 
	double x,y; // Real and imaginary parts

public:
	cpxnumber() : x(0), y(0) {}
	cpxnumber(double x, double y) : x(x), y(y) {}
	cpxnumber(double phase) {
		x = cos(phase);
		y = sin(phase);
	}
	cpxnumber(const cpxnumber &a) : x(a.x), y(a.y) {}
	~cpxnumber () {};

	void set_x(double x);
	void set_y(double y);
	double get_x();
	double get_y();

	void mult_cpx(cpxnumber a);
	void mult_real(double a);
	void add(cpxnumber a);
	void substract(cpxnumber a);
	void print();
	void add(cpxnumber a,cpxnumber b); 
	void substract(cpxnumber a,cpxnumber b); 
	double norm();

};
#endif