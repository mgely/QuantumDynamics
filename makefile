make:
	g++ -o program main.cpp header.h simulation.cpp simulation.h -lm -lfftw3
	python setup.py
	./program
	python plot.py


