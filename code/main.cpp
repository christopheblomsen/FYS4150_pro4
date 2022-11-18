// Build command
// g++ -I include src/* main.cpp -o main -larmadillo
//
// Run command
// ./main output_filename input_name.csv
#include "2D_Ising.hpp"
#include <iostream>
#include <sstream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <vector>


int main(int argc, char** argv){
    std::string mode;
    const char* filename = argv[1];
    bool order;
    arma::mat input;

    if (argc < 2){
        std::cout << "This program should have 2 command line parameters" << std::endl;
        std::cout << "argv[1] is filename where data is saved" << std::endl;
        std::cout << "argv[2] is the input parameter file" << std::endl;
        exit(EXIT_FAILURE);
    }

    std::string input_file = argv[2];
    input.load(arma::csv_name(input_file));
    input.shed_row(0);

    if (input(7) == 0){
        order = false;
    }
    else{
        order = true;
    }

    int L = input(0);
    double J = input(1);
    double N = input(2);
    int cycle = input(3);
    int burn = input(4);
    std::vector<double> temperature = {input(5), input(6)};

    for (double temp : temperature){
        std::string temp_string = std::to_string(temp);
        Ising ise(L, temp, order);

        // SEPARATE CHECK
        // arma::mat lattice = ise.lattice_init(L);
        // std::cout << lattice << std::endl;
        // std::cout << "L is:" << ise.L << std::endl;
        // double E = ise.energy(lattice, J);
        // std::cout << "E_init:" << E << std::endl;

        // double m = ise.magnetization(lattice);
        // std::cout << "m init:" << m << std::endl;

        // ise.Metropolis(lattice, E, m);


        // mcmc check
        arma::mat data = arma::mat(3, (cycle - burn));
        data = ise.mcmc(burn, cycle, data);
        std::cout << "data: " << std::endl << data << std::endl;
        data.save(filename+temp_string+".bin");
    }
    return 0;
}
