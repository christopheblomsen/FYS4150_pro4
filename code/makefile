all: erase compile run plot

erase: 
	rm -f *bin main

compile: 
	g++ -O3 main.cpp -I include src/2D_Ising.cpp -o main -larmadillo -fopenmp 

run:
	./main

plot: 
	python3 plot.py