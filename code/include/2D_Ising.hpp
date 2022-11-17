#include <armadillo>

struct Ising{
    // Temperature is in J/kb so kb should be 1, not 1.380658e-23
    double k{1.};
    int N;
    int L;
    double T;
    double beta;
    Ising(int N, double T);
    arma::mat lattice_init(int N);

    // reset for each mcmc
    // What goes in here?
    void reset();

    // Function for periodic boundary
    int index(int i);

    double energy(arma::mat lattice, int J=1, bool by_spin = false);

    double magnetization(arma::mat lattice, bool by_spin = false);

    double Cv(double E_avg, int N, double T);

    double Chi(arma::vec M, double M_avg, int N, double T);

    arma::mat mcmc(int N_burn, int i, arma::mat data);

    void Metropolis(arma::mat &lattice, double &E_sys, double &M_sys, double &E_avg, double &M_avg);
};
