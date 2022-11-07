#include "2D_Ising.hpp"
#include "omp.h"
#include <iostream>
#include <armadillo>


Ising::Ising(int L, int N_cycle){
    L = L;
    N_cycle = N_cycle;
    N = L*L;
}

void Ising::reset(double* M_tot, double* M_tot2, double* M_abs){
    S = make_matrix(&M_sys);
    (*M_tot)  = 0;
    (*M_tot2) = 0;
    (*M_abs)  = 0;

// create an NxN lattice with random spins from {-1,1}
arma::mat Ising::lattice_init(int N){
    
    //Initialize all spin values
    arma::vec spin = {-1, 1};
    
    // Random seed generated to make sure we never start with the same lattice
    arma::arma_rng::set_seed_random();

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

int Ising::index(int i){
    return (i + L)%L;
}


double Ising::energy(arma::mat lattice, int J, bool by_spin){
    // by_spin if we want to normalize it by spin set to false by default
    int N = lattice.n_rows;
    double E = 0.;
    
    for (int i=0; i < N; i++){
        for (int j=0; j < N; j++){
            E += lattice(i, j) * lattice(i, (j+1)%N) + lattice((i+1)%N, j) * lattice(i, j);
            // std::cout << "J+1: " << (j+1)%N << "  i+1: " << (i+1)%N << std::endl;
            // std::cout << E << std::endl;
        }
    }
    if (by_spin){
        return E *= -J/N;
    }
    else{
        return E *= -J;
    }
}


double Ising::magnetization(arma::mat lattice, bool by_spin){
    // Calculate the magnetization |M| or |m| depending on by_spin, by default
    // is |M|

    if (by_spin){
        int N = lattice.n_cols;
        double M = arma::accu(lattice) / N;
    }

    double M = arma::accu(lattice);

    return M;
}

double Ising::Cv(arma::vec E, int N, double T){
    
    arma::vec E2 = arma::pow(E, 2);
    double average_E = arma::mean(E);
    double average_E2 = arma::mean(E2);

    double cv = (average_E2 - average_E * average_E2) / (N * T * T);
    return cv;

}


double Ising::Chi(arma::vec M, int N, double T){
    
    // make the 2 version of M needed
    arma::vec M2 = arma::pow(M, 2);
    arma::vec absM = arma::abs(M);

    // take the average
    double average_M2 = arma::mean(M2);
    double average_absM = arma::mean(absM);


    //Now compute the susceptibility
    double xi = (average_M2 - average_absM * average_absM) / (N * T);

    return xi;
}

void Ising::mcmc(int N_burn, int i, arma::vec* Cv_vec){
    int N = L*L;
    reset(&M_tot, &M_tot2, &M_abs);


}

void Ising::Metropolis(arma::mat &lattice, double* E_sys, double* M_sys){
    for (int n = 0; n < N; n++){
        int i = arma::randi(arma::distri_param(0, L-1));
        int j = arma::randi(arma::distri_param(0, L-1));

        int dE = 2*lattice(index(i), index(j))*
            (lattice(index(i+1), index(j))
             + lattice(index(i-1), index(j))
             + lattice(index(i), index(j+1))
            + lattice(index(i), index(j-1)));

        if (dE <= 0){
            lattice(i, j) *= -1;
            int dM = 2*lattice(i, j);
            (*E_sys) += dE;
            (*M_sys) += dM;
        }
        else{
            lattice(i, j) *= -1;
            int dM = 2*lattice(i, j);
            (*E_sys) += dE;
            (*M_sys) += dM;
        }
    }
}
