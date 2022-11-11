#include "2D_Ising.hpp"
#include "omp.h"
#include <iostream>
#include <armadillo>


Ising::Ising(int N, double T){
    L = L;
    N_cycle = N_cycle;
    N = L*L;
    beta = 1/(k*T);
}

// No idea
void Ising::reset(){
    std::cout << "No idea" << std::endl;
}

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
    // if (by_spin){
    //     return E *= -J/N;
    // }
    // else{
    //     return E *= -J;
    // }
    return E;
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

double Ising::Cv(double E_avg, int N, double T){
    // Avg E and avg E ^2

//     // arma::vec E2 = arma::pow(E, 2);
//     // double average_E = arma::mean(E);
//     // double average_E2 = arma::mean(E2);


    double cv = (E_avg - E_avg * E_avg) / (N * T * T);
    return cv;

}


double Ising::Chi(double M_avg, int N, double T){
    double xi = (M_avg - M_avg * M_avg) / (N * T);
    return xi;
}

void Ising::mcmc(int N_burn, int i, arma::vec Cv_vec){
    arma::mat lattice = lattice_init(N);
    double E_sys = energy(lattice);
    double M_sys = magnetization(lattice);

    double E_sum = 0.;
    double M_sum = 0.;

    for (int i=0; i < N_burn; i++){
        Metropolis(lattice, E_sys, M_sys);
        E_sum += E_sys;
        M_sum += std::fabs(M_sys);
    }
    double E_avg = E_sum/N_burn;
    double M_avg = M_sum/N_burn;
}

void Ising::Metropolis(arma::mat &lattice, double &E_sys, double &M_sys){
    for (int n = 0; n < N; n++){

        // Pick a random spin
        int i = arma::randi(arma::distr_param(0, L-1));
        int j = arma::randi(arma::distr_param(0, L-1));

        int dE = 2*lattice(index(i), index(j))*
            (lattice(index(i+1), index(j))
             + lattice(index(i-1), index(j))
             + lattice(index(i), index(j+1))
            + lattice(index(i), index(j-1)));

        if (dE <= 0){
            lattice(i, j) *= -1;
            int dM = 2*lattice(i, j);
            (E_sys) += dE;
            (M_sys) += dM;
        }
        else{
            if (arma::randu() < std::exp(-beta*dE)){
                lattice(i, j) *= -1;
                int dM = 2*lattice(i, j);
                (E_sys) += dE;
                (M_sys) += dM;
            }
        }
    }
}
