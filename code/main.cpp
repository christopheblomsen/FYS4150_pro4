#include "2D_Ising.hpp"
#include <iostream>


int main(){
    
    int N = 5;
    
    arma::mat lattice = lattice_init(N);
    std::cout << lattice << std::endl;

    return 0;
}