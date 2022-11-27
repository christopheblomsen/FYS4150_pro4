#ifndef UTILS_H_
#define UTILS_H_

#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <vector>
#include "omp.h"
#include "chrono"

// runs problem 4
void problem4();

// runs problem 5
void problem5();

// runs problem 7
void problem7(int L=20, int burn=0,
              int cycle=50000, bool order=false);

// runs problem 8
void problem8();

#endif // UTILS_H_
