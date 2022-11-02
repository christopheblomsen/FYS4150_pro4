#include "2D_Ising.hpp"
#include <iostream>


int main(){
    
    int N = 3;
    double J = 1.;
    
    arma::mat lattice = lattice_init(N);
    std::cout << lattice << std::endl;

    double E = energy(lattice, J);
    std::cout << E << std::endl;

    return 0;
}