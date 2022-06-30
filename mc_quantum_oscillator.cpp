//#include "mpi.h"
#include<iostream>
#include<cmath>
#include<vector>
#include<random>
#include<cstdlib>
#include<algorithm>
#include <omp.h>

void Metropolis_Hasting_Serial(std::vector<double> &U, int moves, int N, int seed, double a, double g, double omega, double m, double eps);
/*std::vector<double> Metropolis_Hasting_Parallel(std::vector<double> &U, int movesh, int N, int seed, double a, double g, int pid, int np);
 */

//auxiliary functions
double action_change( std::vector<double> &U, double shift, int N, int t, double g);
double transform_u_to_x(double u, const double M, const double EPS, const double OMEGA);
double x_sq_expected_value(std::vector<double> U, int gap, int jj, int &counter, int &X_sq, int m, int eps, int omega, int N);

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
  std::uniform_real_distribution<double> filling(-100.0,100.0);
  
  //prueba
  
    std::vector <double> U; U.resize(N); //Inicialización de la variable U
  //U.resize(Nt);
    
  //inicialización del problema

  std::fill(U.begin(), U.end(), 0.0);  //cold start
  //std::fill(U.begin(), U.end(), filling(gen));   //hot start
  U[0]=0.0; U[N-1]=0.0; //condiciones iniciales	

  Metropolis_Hasting_Serial(U,  moves, N, rd(), a, g, omega, m, epsilon);

  
  //MPI_setup
  /*int np, pid;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &np);
  MPI_Comm_rank(MPI_COMM_WORLD, &pid);

  
  int Nt=N/np;
  std::vector <double> U; U.resize(N); //Inicialización de la variable U
  //U.resize(Nt);
    
  //inicialización del problema

  std::fill(U.begin(), U.end(), 0.0);  //cold start
  //std::fill(U.begin(), U.end(), filling(gen));   //hot start
  U[0]=0.0; U[N-1]=0.0; //condiciones iniciales	

  Metropolis_Hasting_Serial(U,  moves, N, seed, a, g);
 
  MPI_Finalize();
 
  */
  return 0;
}

/*
std::vector<double> Metropolis_Hasting_Parallel(std::vector<double> &U, int moves, int N, int seed, double a, double g, int pid, int np){

  std::mt19937 gen(seed);
  std::uniform_real_distribution<double> dis(-a,a);
  std::uniform_real_distribution<double> P(0,1);
  int  t=0; //time slice

  //parallel setup

  double Nt= N/np;

  std::vector <double> U_init; U_init.resize(N);
  std::vector <double> U_new; U_new.resize(N);

  U_init=U;
  
  //Variables a calcular
  double ratio=0.0;
  double Delta_S=0.0;  
  
  for(int jj=0; jj<Nt;++jj){
    for(int ii=0; ii<moves;++ii){
      t=N*std::floor(P(gen));
      if(t==0 || t==N-1){
	break;
      }
      else{
	U_new[t]=U_init[t]+dis(gen);
	double alpha=P(gen);
	Delta_S=action_change(U_init, U_new, N, t, g);
	ratio=std::min(1.0, std::exp(-Delta_S));
	if(ratio>=alpha){
	  U[t]=U_new[t];
	}
      }
    }
  }
  return U;
}


*/



void Metropolis_Hasting_Serial(std::vector<double> &U, int Nsweeps, int N, int seed, double a, double g, double omega, double m, double eps){

  std::mt19937 gen(seed);
  std::uniform_real_distribution<double> dis(-a,a);
  std::uniform_real_distribution<double> P(0.0,1.0);
  int  t=0; //time slice
  int gap=10;
  double  counter=0.0;
  double X_sq=0.0;

  
  //Variables a calcular
  double ratio=0.0;
  double Delta_S=0.0;  
  double alpha=0;
  double shift=0.0;
  
 
  
  for(int jj=1; jj<=Nsweeps;jj++){
    for(t=1; t<N-1;t++){         //medición tomando cada uno de los t's
      //t=std::floor((N-1)*P(gen));     //medicion tomando un t aleatorio, se probo ambas opciones,  arrojando resultados que se acercan al valor esperado en ambos, pero con mas exactitud para el caso de arriba
      // if(t!=0 && t!=N-1){
	shift=dis(gen);
	alpha=P(gen); 
	Delta_S=action_change(U, shift, N, t, g); //cambio de la acción
	ratio=std::min(1.0, std::exp(-Delta_S));
	if(alpha<=ratio){
	  U[t] += shift;
	}	
	//} 
    }
    /* if (jj%gap==0){
      ++counter;
      for (int k=0; k<N; ++k){
	X_sq += std::pow( transform_u_to_x(U[k], m,eps, omega), 2);
      }
      //double mean_sq= X_sq/(N*counter);
	//std::cout<<jj <<"  "<<mean_sq<<std::endl;
    }
  }
  double mean_sq= X_sq/(N*counter);*/
    std::cout<<"The expected value of x² is "<< x_sq_expected_value(U, gap, jj, counter, X_sq, m, eps, omega, N)<<std::endl;
  }
}


//auxiliary functions

double action_change( std::vector<double> &U, double shift, int N, int t, double g){
  
  double ds = 0.0;
   ds= shift*(shift+2.0*U[t]-g*(U[t-1]+U[t+1]));
  return ds;
}


double transform_u_to_x(double u, const double M, const double EPS, const double OMEGA){
  double k = std::pow(EPS * OMEGA *0.5, 2);
  double x = std::sqrt(EPS/(M*(k+1)))*u;
  return x;
}

double x_sq_expected_value(std::vector<double> U, int gap, int jj, int &counter, int &X_sq, int m, int eps, int omega, int N){
  //no he podido modulariar bien esto :c
  if (jj%gap==0){
    ++counter;
    for (int k=0; k<N; ++k){
      X_sq += std::pow( transform_u_to_x(U[k], m,eps, omega), 2);
    }
  }
  double mean_sq= X_sq/(N*counter);

  return mean_sq;
}
