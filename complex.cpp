#include "header.h"

void cpxnumber::set_x(double in){x=in;}
void cpxnumber::set_y(double in){y=in;}
double cpxnumber::get_x(){return x;}
double cpxnumber::get_y(){return y;}

void cpxnumber::mult_cpx(cpxnumber a){
	
	double x_temp = x;
	x = a.x*x - a.y*y;
	y = a.x*y + a.y*x_temp;
}

void cpxnumber::mult_real(double a){
	x = a*x;
	y = a*y;
}

void cpxnumber::add(cpxnumber a){
	x = a.x + x;
	y = a.y + y;
}

void cpxnumber::substract(cpxnumber a){
	x = x - a.x ;
	y = y - a.y;
}

void cpxnumber::print(){
	std::cout<<"Re part: " << x <<std::endl;
	std::cout<<"Im part: " << y <<std::endl;
}

void cpxnumber::add(cpxnumber a,cpxnumber b){
	x = a.x + b.x;
	y = a.y + b.y;
}

void cpxnumber::substract(cpxnumber a,cpxnumber b){
	x = a.x - b.x;
	y = a.y - b.y;
}

double cpxnumber::norm(){
	return x*x+y*y;
}