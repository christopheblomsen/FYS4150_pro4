#include "2D_Ising.hpp"
#include <iostream>
#include <armadillo>
#include <random>


// create an NxN lattice with random spins from {-1,1}
arma::mat lattice_init(int N){
    
    //Initialize all spin values
    arma::vec spin = {-1, 1};
    
    // Random seed generated to make sure we never start with the same lattice
    srand(time(NULL));

    arma::mat lattice = arma::mat(N, N);
    // Apparently setting all initial random values at once is faster
    arma::mat random_value = arma::randi<arma::mat>(N, N, arma::distr_param(0, 1));
    std::cout << random_value << std::endl;


    for (int i=0; i < N; i++){
        for (int j=0; j < N; j++){

            int index = random_value(i, j);
            lattice(i, j) = spin(index);
        }

    }

    return lattice;
}


double energy(arma::mat lattice, int J){

    int N = lattice.n_rows;
    double E = 0.;
    
    for (int i=0; i < N; i++){
        for (int j=0; j < N; j++){
            E += lattice(i, j) * lattice(i, (j+1)%N) + lattice((i+1)%N, j) * lattice(i, j);
            // std::cout << "J+1: " << (j+1)%N << "  i+1: " << (i+1)%N << std::endl;
            // std::cout << E << std::endl;
        }
    }
    return E *= -J;
}