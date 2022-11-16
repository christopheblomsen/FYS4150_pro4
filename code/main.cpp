#include "2D_Ising.hpp"
#include <iostream>


int main(){
    
    int L = 2;
    double J = 1.;
    double N = 4;
    double temperature = 1.;
    int cycle = 2;
    int burn = 0;

    Ising ise(L, temperature);

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
    arma::mat data = arma::mat(2, (cycle - burn));
    data = ise.mcmc(burn, cycle, data);
    std::cout << "data: " << data << std::endl;     
    return 0;
}
