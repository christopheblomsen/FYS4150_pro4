// Build command
// g++ -I include src/* main.cpp -o main -larmadillo -fopenmp
//
// Run command
// ./main N
//
// Where N is the problem number
#include "2D_Ising.hpp"
#include "utils.h"

int main(int argc, char** argv){
    int prob;
    if (argc == 1){
        prob = 0;
    }
    else{
        prob = std::atoi(argv[1]);
    }

    switch(prob){
        case 4:
            problem4();
            break;
        case 5:
            problem5();
            break;
        case 7:
            problem7();
            break;
        case 8:
            problem8();
            break;
        default:
            std::cout << "No problem argument given" << std::endl;
            break;
    }
    return 0;
}

 // for (double temp : temperature){
    //     std::string temp_string = std::to_string(temp);
    //     Ising ise(L, temp, order);

    //     // SEPARATE CHECK
    //     // arma::mat lattice = ise.lattice_init(L);
    //     // std::cout << lattice << std::endl;
    //     // std::cout << "L is:" << ise.L << std::endl;
    //     // double E = ise.energy(lattice, J);
    //     // std::cout << "E_init:" << E << std::endl;

    //     // double m = ise.magnetization(lattice);
    //     // std::cout << "m init:" << m << std::endl;

    //     // ise.Metropolis(lattice, E, m);


    //     // mcmc check
    //     // arma::mat data = arma::mat(3, (cycle - burn));
    //     // data = ise.mcmc(burn, cycle, data);
    //     // std::cout << "data: " << std::endl << data << std::endl;
    //     // data.save(filename+temp_string+".bin");
    // }
