#include "utils.h"
#include "2D_Ising.hpp"
#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <vector>
#include <armadillo>

// Runs problem 4
void problem4(){
    std::cout << "Running problem 4" << std::endl;
    double T_1 = 1.;
    int L_1 = 2;

    Ising ising1(L_1, T_1);
    arma::vec cycles_1 = {500, 1000 ,2000, 3000, 4000, 5000, 6000, 7000,
    8000, 9000, 10000, 20000, 30000, 40000, 50000, 60000, 70000, 80000, 90000, 100000};

    for (int cycle : cycles_1){
        arma::mat data = arma::mat(3, cycle);
        data = ising1.mcmc(0, cycle, data);
        data.save("L2_T1_"+std::to_string(cycle)+"_cycles.bin");
    }

}

// Runs problem 5
void problem5(){
    std::cout << "Running problem 5" << std::endl;
    // Set up the parameters for the test
    int L_2 = 20;
    double T_2 = 1.;
    double T_21 = 2.4;

    // Create the environnement for both temperature
    Ising ising2(L_2, T_2);
    Ising ising2_1(L_2, T_21);
    Ising ising2_ord(L_2, T_2, true);
    Ising ising2_1_ord(L_2, T_21, true);

    arma::vec cycles_1 = {500, 1000 ,2000, 3000, 4000, 5000, 6000, 7000,
    8000, 9000, 10000, 20000, 30000, 40000, 50000, 60000, 70000, 80000, 90000, 100000};

    // Iterate over every number of cycle
    for (int cycle : cycles_1){

        // Create four data different for both temperature and start
        arma::mat data1 = arma::mat(3, cycle);
        arma::mat data2 = arma::mat(3, cycle);
        arma::mat data1_ord = arma::mat(3, cycle);
        arma::mat data2_ord = arma::mat(3, cycle);

        data1 = ising2.mcmc(0, cycle, data1);
        data2 = ising2_1.mcmc(0, cycle, data2);
        data1_ord = ising2_ord.mcmc(0, cycle, data1);
        data2_ord = ising2_1_ord.mcmc(0, cycle, data2_ord);

        data1.save("L20_random_T1_"+std::to_string(cycle)+"_cycles.bin");
        data2.save("L20_random_T24_"+std::to_string(cycle)+"_cycles.bin");
        data1_ord.save("L20_ord_T1_"+std::to_string(cycle)+"_cycles.bin");
        data2_ord.save("L20_ord_T24_"+std::to_string(cycle)+"_cycles.bin");
    }
}

// Runs problem 7
void problem7(int L, int burn, int cycle, bool order){
    std::cout << "Running problem 7" << std::endl;
    // MCMC parallel check
    arma::vec temp = arma::linspace(1., 2.4, 50);
    temp.save("temperature.bin");

    auto start1 = std::chrono::high_resolution_clock::now();
    // parallize the for loop so over every temperature
    #pragma omp parallel for
    for (int t=0; t < temp.n_elem; t++){


        // create the data

        arma::mat data = arma::mat(3, (cycle - burn));
        double T = temp[t];

        // create the environment and run the simulation
        Ising ise(L, T, order);
        data = ise.mcmc(burn, cycle, data);
        // std::cout << T <<std::endl;
        // data.save("test_para"+std::to_string(t)+".bin");
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
        data = ise.mcmc(burn, cycle, data);
        // std::cout << T <<std::endl;
        // data.save("test_serial_"+std::to_string(t)+".bin");
        // std::cout << "temperature: " << T << std::endl;
    }
    auto finish2 = std::chrono::high_resolution_clock::now();

    std::cout << "Time for serial: " << std::chrono::duration_cast<std::chrono::seconds>(finish2-start2).count() << " s" << std::endl;

}

// Runs problem 8
void problem8(){
    std::cout << "Running problem 8" << std::endl;
    arma::vec L_list = {40, 60, 80, 100};
    arma::vec temperature8 = arma::linspace(2.1, 2.4, 10);
    int cycle8 = 100000;
    int burn = 50000;

    temperature8.save("temp_prob8.bin");

    #pragma omp parallel for collapse(2)
    for (int t=0; t < temperature8.n_elem; t++){
        for (int L : L_list){
            double T8 = temperature8[t];
            // We create the data
            arma::mat data = arma::mat(3, (cycle8 - burn));
            // We create the environment and run the simulation
            Ising ise(L, T8);
            data = ise.mcmc(burn, cycle8, data);
            data.save("L"+std::to_string(L)+"_prob8_"+std::to_string(t)+".bin");
        }
    }
}
