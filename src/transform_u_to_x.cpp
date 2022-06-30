#include "transform_u_to_x.hpp"

double transform_u_to_x(double u, const double M, const double EPS, const double OMEGA){
  double k = std::pow(EPS * OMEGA *0.5, 2);
  double x = std::sqrt(EPS/(M*(k+1)))*u;
  return x;
}

