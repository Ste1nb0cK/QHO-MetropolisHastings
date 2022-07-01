#include "mpi.h"
#include "action_change.hpp"
#include "transform_u_to_x.hpp"
#include "fill.hpp"
#include "Metropolis_Hasting_Parallel.hpp"
#include "Metropolis_Hasting_Serial.hpp"
#include "x_sq_Metropolis_Hasting_Serial.hpp"
#include<iostream>
#include<cmath>
#include<vector>
#include<random>
#include<cstdlib>
#include<algorithm>

int main(int argc, char **argv){

  //Variables del problema
  const int  N= std::atoi(argv[1]); //Discretización del tiempo
  const int moves= std::atoi(argv[2]); //Cantidad de movimientos en el sampleo
  const double epsilon=0.25;
  const double omega=1;
  const double k=std::pow(epsilon*omega*0.5,2);
  const double g=(1-k)/(1+k);
  const double m=std::atoi(argv[3]);

		
  //Distribuciones necesarias
  std::random_device rd;
  const double a=0.5;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<double> dis(-a,a);
  std::uniform_real_distribution<double> P(0,1);
  
  //MPI_setup
  int np, pid;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &np);
  MPI_Comm_rank(MPI_COMM_WORLD, &pid);

  int Nt=N/np;
  std::vector <double> U;  //Inicialización de la variable U
  //U.resize(N); //size for serial version
  U.resize(Nt); //size for parallel version
    


  Metropolis_Hasting_Parallel(U, moves,  N, rd(), a,  g,  omega,  m,  epsilon,  pid,  np);
  
  MPI_Finalize();
 
 
  return 0;
}


