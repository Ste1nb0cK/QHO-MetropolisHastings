#pragma once
#include "mpi.h"
#include<iostream>
#include<cmath>
#include<vector>
#include<random>
#include<cstdlib>
#include<algorithm>
#include "action_change.hpp"
#include "transform_u_to_x.hpp"

void Metropolis_Hasting_Serial(std::vector<double> &U, int Nsweeps, int N, int seed, double a, double g, double omega, double m, double eps);
