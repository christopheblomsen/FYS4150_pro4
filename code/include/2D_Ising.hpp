#include <armadillo>

struct Ising{
    // Temperature is in J/kb so kb should be 1, not 1.380658e-23
    double k{1.};
    int N;
    int L;
    double T;
    double beta;
    bool order;
    Ising(int length, double temperature, bool order_=false);
    arma::mat lattice_init(int N);

    // Method for ordered lattice
    void ordered_lattice(int N, int spin, arma::mat &lattice);

    void unordered_lattice(int N, arma::vec spin, arma::mat &lattice);

    // Function for periodic boundary
    int index(int i);

    double energy(arma::mat &lattice, int J=1, bool by_spin = false);

    double magnetization(arma::mat &lattice, bool by_spin = false);

    double Cv(arma::vec E, int N, double T);

    double Chi(arma::vec M, int N, double T);

    arma::mat mcmc(int N_burn, int i, arma::mat data);

    void Metropolis(arma::mat &lattice, double &E_sys, double &M_sys);
};
