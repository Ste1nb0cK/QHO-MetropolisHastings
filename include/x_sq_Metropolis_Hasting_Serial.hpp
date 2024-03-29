#pragma once

#include<iostream>
#include<cmath>
#include<vector>
#include<random>
#include<cstdlib>
#include<algorithm>

double x_sq_Metropolis_Hasting_Serial(std::vector<double> &U, int moves, int N, int seed, double a, double g, double omega, double m, double eps);
