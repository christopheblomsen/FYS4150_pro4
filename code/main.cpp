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
#include "omp.h"
#include "chrono"


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

    // MCMC parallel check
    arma::vec temp = arma::linspace(temperature[0], temperature[1], 100);

    auto start1 = std::chrono::high_resolution_clock::now();
    // parallize the for loop so over every temperature 
    #pragma omp parallel for
    for (int t=0; t < temp.n_elem; t++){
        // create the data 
        arma::mat data = arma::mat(3, (cycle - burn));
        double T = temp[t];

        // create the environment and run the simulation
        Ising ise(L, T, order);
        data = ise.mcmc(cycle, burn, data);
        // data.save("test_2"+std::to_string(T)+".bin");
    }


    auto finish1 = std::chrono::high_resolution_clock::now();


    std::cout << "Time for parallel: " << std::chrono::duration_cast<std::chrono::seconds>(finish1-start1).count() << " s" << std::endl;

    auto start2 = std::chrono::high_resolution_clock::now();
    //Without parallelization
    std::cout << "elem: " << temp.n_elem <<std::endl;

    for (int t=0; t < temp.n_elem; t++){
        // create the data 
        arma::mat data = arma::mat(3, (cycle - burn));
        double T = temp[t];

        // create the environment and run the simulation
        Ising ise(L, T, order);
        data = ise.mcmc(cycle, burn, data);
        // data.save("test_serial"+std::to_string(T)+".bin");
        // std::cout << "temperature: " << T << std::endl;
    }
    auto finish2 = std::chrono::high_resolution_clock::now();

    std::cout << "Time for serial: " << std::chrono::duration_cast<std::chrono::seconds>(finish2-start2).count() << " s" << std::endl;


    // for (double temp : temperature){
    //     std::string temp_string = std::to_string(temp);
    //     Ising ise(L, temp, order);

    //     // SEPARATE CHECK
    //     // arma::mat lattice = ise.lattice_init(L);
    //     // std::cout << lattice << std::endl;
    //     // std::cout << "L is:" << ise.L << std::endl;
    //     // double E = ise.energy(lattice, J);
    //     // std::cout << "E_init:" << E << std::endl;

    //     // double m = ise.magnetization(lattice);
    //     // std::cout << "m init:" << m << std::endl;

    //     // ise.Metropolis(lattice, E, m);


    //     // mcmc check
    //     // arma::mat data = arma::mat(3, (cycle - burn));
    //     // data = ise.mcmc(burn, cycle, data);
    //     // std::cout << "data: " << std::endl << data << std::endl;
    //     // data.save(filename+temp_string+".bin");
    
    // }
    return 0;
}
