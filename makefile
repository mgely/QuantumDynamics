make:
	g++ -o program main.cpp header.h simulation.cpp simulation.h initialization.cpp -lm -lfftw3
	./program
	rm plot.dat gnusettings2d.txt program

