all: erase compile run plot

erase: 
	rm -f *bin main

compile: 
	g++ main.cpp -I include src/2D_Ising.cpp -o main -larmadillo 

run:
	./main

plot: 
	python3 plot.py