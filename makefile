make:
	g++ -o program main.cpp header.h quantum.cpp quantum.h FFT.cpp FFT.h complex.h complex.cpp
	./program
	python data.py


