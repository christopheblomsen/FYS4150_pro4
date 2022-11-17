#include "2D_Ising.hpp"
#include "omp.h"
#include <iostream>
#include <armadillo>


Ising::Ising(int length, double temperature, bool order_){
    L = length;
    N = L*L;
    T = temperature;
    beta = 1/(k*T);
    order = order_;

    // Probably could add the constant for our system 
    // when we calculate Cv and Chi 
}

// No idea,  maybe just recreate a new lattice? 
void Ising::reset(){
    std::cout << "No idea" << std::endl;
}

void Ising::ordered_lattice(int N, int spin, arma::mat &lattice){
    for (int i=0; i < N; i++){
        for (int j=0; j < N; j++){
            lattice(i, j) = spin;
        }
    }
}

void Ising::unordered_lattice(int N, arma::vec spin, arma::mat &lattice){
    // Random seed generated to make sure we never start with the same lattice
    arma::arma_rng::set_seed_random();
    // Apparently setting intall initial random values at once is faster
    arma::mat random_value = arma::randi<arma::mat>(N, N, arma::distr_param(0, 1));
    // std::cout << random_value << std::endl;

    for (int i=0; i < N; i++){
        for (int j=0; j < N; j++){
            int index = random_value(i, j);
            lattice(i, j) = spin(index);
        }
    }
}
// create an NxN lattice with random spins from {-1,1}
arma::mat Ising::lattice_init(int N){
    
    //Initialize all spin values
    arma::vec spin = {-1, 1};
    arma::mat lattice = arma::mat(N, N);
    
    if (order){
        ordered_lattice(N, 1, lattice);
    }
    else{
        unordered_lattice(N, spin, lattice);
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

    // For the 2x2 case we count twice the energy with this method
    if (N == 2){
        E = E /2.;
    }

    if (J == 1){
        return E;
    }
    else{
        return E*J;
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


    double cv = (average_E2 - average_E * average_E) / (N * T * T);
    return cv;

}


double Ising::Chi(arma::vec M, int N, double T){
    // Remove M_avg signature
    arma::vec M_ = arma::pow(M, 2);
    double M2 = arma::mean(M_);
    double M_avg = arma::mean(M);

    double xi = (M2 - M_avg * M_avg) / (N * T);
    return xi;
}


// We run several cycles of Metropolis  
arma::mat Ising::mcmc(int N_burn, int i, arma::mat data){
    // N_burn is the number of cycle we want to throw away 
    // i is the number of cycle we want to do 
    // data is matrix to get the energy and the magnetisation
    // the values of interest are calculated outside to help 
    // better recognize where a potential bug comes from. 

    // Initialize the system 
    arma::mat lattice = lattice_init(L);
    double E_sys = energy(lattice);
    double M_sys = magnetization(lattice);

    // // Save the initial condition, i.e energy and magnetization
    // data(0, 0) = E_sys;
    // data(1, 0) = M_sys;

    // To keep track in our data matrix index
    int index = 0;

    // We run the i cycles 
    for (int k=1; k <= i; k++){
        
        // We reset the averaged value to the initial system 
        // after every cycle
        double E_avg = E_sys;
        double M_avg = M_sys;

        Metropolis(lattice, E_sys, M_sys, E_avg, M_avg);

        // When we pass the burning point we start registering the values
        if (k > N_burn){
            data(0, index) = E_sys;
            data(1, index) = std::fabs(M_sys);
            data(2, index) = k;

            index += 1; 
        }
    }
    
    return data;
    // double E_sum = 0.;
    // double M_sum = 0.;

    // for (int i=0; i < N_burn; i++){
    //     Metropolis(lattice, E_sys, M_sys);
    //     E_sum += E_sys;
    //     M_sum += std::fabs(M_sys);
    // }
    // double E_avg = E_sum/N_burn;
    // double M_avg = M_sum/N_burn;
}


// This is to run one iteration of the Metropolis algorithm
// modify the lattice and update the energy and magnetization (only the current step)
void Ising::Metropolis(arma::mat &lattice, double &E_sys, double &M_sys,
                       double &E_avg, double &M_avg){
    
    //  One cycle correspond to N spin flip attempt. 
    for (int n = 0; n < N; n++){
        std::cout << "Iteration: " << n+1 << std::endl;

        // Pick a random spin
        int i = arma::randi(arma::distr_param(0, L-1));
        int j = arma::randi(arma::distr_param(0, L-1));
        std::cout << "i: " << i << std::endl;
        std::cout << "j: " << j << std::endl;


        // to Calculate dE we first compute the energy around 
        // the proposed spin and the energy around if flipped
        // the soustraction gives us the deltaE 
        int E_noflip = lattice(index(i), index(j))*
            (lattice(index(i+1), index(j))
             + lattice(index(i-1), index(j))
             + lattice(index(i), index(j+1))
            + lattice(index(i), index(j-1)));

        int E_flip = -1 * E_noflip;
        int dE = E_flip - E_noflip;

        // If we have a 2x2 lattice we 
        if (L == 2){
            dE /= 2.;
        }

        // Should always be -8, -4, 0, 4, 8 
        std::cout << "dE is: " << dE << std::endl;
        std::cout << lattice << std::endl;
        std::cout << "entering step verification" << std::endl;
        std::cout << "" << std::endl;

        // Verification of the new state. 
        if (dE <= 0){
            lattice(i, j) *= -1;
            int dM = 2*lattice(i, j);
            (E_sys) += dE;
            (M_sys) += dM;

            std::cout << "Change done (<=0)" << std::endl;
            std::cout << "New energy is: " << E_sys << std::endl;
            std::cout << "New M is: " << M_sys << std::endl; 
            std::cout << lattice << std::endl;
            std::cout << "" << std::endl;
        
        }
        else{
            double u = arma::randu();
            double p = std::exp(-beta * dE);

            std::cout << "u is: " << u << std::endl;
            std::cout << "p is: " << p << std::endl;

            if (u < p){
                lattice(i, j) *= -1;
                int dM = 2*lattice(i, j);
                (E_sys) += dE;
                (M_sys) += dM;

                std::cout << "Change done" << std::endl;
                std::cout << "New energy is: " << E_sys << std::endl;
                std::cout << "New M is: " << M_sys << std::endl; 
                std::cout << lattice << std::endl;
                std::cout << "" << std::endl;
            }
        std::cout << "No change" << std::endl;
        std::cout << "energy is: " << E_sys << std::endl;
        std::cout << "M is: " << M_sys << std::endl;
        std::cout << lattice << std::endl;
        std::cout << "" << std::endl;

        }

        // Even said useless
        // // update the average
        // E_avg += E_sys;
        // M_avg += M_sys;
    }

    // Take the average
    // E_avg = E_avg / (N + 1);
    // M_avg = M_avg / (N + 1);
}
// Ordered and random in init, plotting in python
//
// prob 6, sample the mcmc. for T=1, only 3(4) different energy
// Not to wide, not wider then something
//
//
// prob 8, many different temperature, energy as function of temp
