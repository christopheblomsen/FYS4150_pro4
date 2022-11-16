#include "2D_Ising.hpp"
#include <iostream>


int main(){
    
    int L = 3;
    double J = 1.;
    double N = 4;

    Ising ise(L, 1);
    arma::mat lattice = ise.lattice_init(L);
    std::cout << lattice << std::endl;
    std::cout << "L is:" << ise.L << std::endl;
    double E = ise.energy(lattice, J);
    std::cout << "E_init:" << E << std::endl;

    double m = ise.magnetization(lattice);
    std::cout << "m init:" << m << std::endl;

    ise.Metropolis(lattice, E, m);

    return 0;
}
