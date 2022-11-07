#include <armadillo>

struct Ising{
    int N;
    int N_cycle;
    int L;
    Ising(int N_cycle);
    arma::mat lattice_init(int L, int N_cycle);

    // reset for each mcmc
    void reset(double* M_tot, double* M_tot2, double* M_abs);

    // Function for periodic boundary
    int index(int i);

    double energy(arma::mat lattice, int J, bool by_spin = false);

    double magnetization(arma::mat lattice, bool by_spin = false);

    double Cv(arma::vec E, int N, double T);

    double Chi(arma::vec M, int N, double T);

    void mcmc(int N_burn, int i, arma::vec* Cv_vec);

    void Metropolis(int n_spins, long& idum, int **spin_matrix,
                    double& E, double& M, double *w);
};
