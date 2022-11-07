#include "2D_Ising.hpp"
#include <iostream>


int main(){
    
    int N = 3;
    double J = 1.;
    Ising ise(N);
    arma::mat lattice = ise.lattice_init(N);
    std::cout << lattice << std::endl;

    double E = ise.energy(lattice, J);
    std::cout << E << std::endl;

    double m = ise.magnetization(lattice);
    std::cout << m << std::endl;

    return 0;
}
