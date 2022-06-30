#include "fill.hpp"

void fill (std::vector<double> &U, int N, bool cold){

  std::random_device rd;
  std::mt19937 gen(rd());
  
  if(cold==true){
    std::fill(U.begin(), U.end(), 0.0);  //cold start
  } else{
    for(int ii=0;ii<N;++ii){
      U[ii]= rd();
    }//hot start
  }
  
  U[0]=0.0; U[N-1]=0.0; //condiciones de frontera	
}
