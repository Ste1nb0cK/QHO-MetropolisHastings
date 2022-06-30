#pragma once
#include "mpi.h"
#include<iostream>
#include<cmath>
#include<vector>
#include<random>
#include<cstdlib>
#include<algorithm>

void Metropolis_Hasting_Parallel(std::vector<double> &U, int Nsweeps, int N, int seed, double a, double g, double omega, double m, double eps, int pid, int np);
