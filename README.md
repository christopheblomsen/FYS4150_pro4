# FYS4150 project 4

This repository is for [project 4](https://anderkve.github.io/FYS3150/book/projects/project4.html) in FYS4150. Which takes a numerical approach on the Ising model. The [report](report/report.pdf) and [.tex](report/report.tex) file can be found under the [report](report/) directory. 

The code is found in the [code](code/) directory. And is divided into the [header](code/include/) directory named `include` and [source](code/src/) directory named `src`. To build the code

```Shell
/code $ g++ -I include src/* main.cpp -o main -larmadillo -fopenmp
```

The run command

```Shell
code $ ./main N
```

Where `N` is an integer number representing the problem we want to run. And can take the numbers `4, 5, 7, 8`. 
The corresponding problems can we found in the [utils.cpp](code/utils.cpp) file.
