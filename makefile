make:
	g++ -o program main.cpp header.h simulation.cpp simulation.h initialization.cpp -lm -lfftw3
	# python setup.py
	./program

