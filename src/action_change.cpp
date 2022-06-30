#include "action_change.hpp"

double action_change( std::vector<double> &U, double shift, int N, int t, double g){
  
  double ds = 0.0;
   ds= shift*(shift+2.0*U[t]-g*(U[t-1]+U[t+1]));
  return ds;
}
