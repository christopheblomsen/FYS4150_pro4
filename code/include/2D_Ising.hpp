#include <armadillo>

arma::mat lattice_init(int N);

double energy(arma::mat lattice, int J, bool by_spin = false);

double magnetization(arma::mat lattice, bool by_spin = false);

double Cv(arma::vec E, int N, double T);

double Chi(arma::vec M, int N, double T);