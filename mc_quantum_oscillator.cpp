//#include "mpi.h"
#include<iostream>
#include<cmath>
#include<vector>
#include<random>
#include<cstdlib>
#include<algorithm>
#include <omp.h>


void Metropolis_Hasting_Serial(std::vector<double> &U, int moves, int N, int seed, double a, double g, double omega, double m, double eps);
void Metropolis_Hasting_Serial_2(std::vector<double> &U, int moves, int N, int seed, double a, double g, double omega, double m, double eps);
/*std::vector<double> Metropolis_Hasting_Parallel(std::vector<double> &U, int moves, int N, int seed, double a, double g, int pid, int np);
*/double action_change(std::vector<double> U_init, std::vector<double> U_new, int N, int t, double g);
double transform_u_to_x(double u, const double M, const double EPS, const double OMEGA);
double action_change_2( std::vector<double> U_init, std::vector<double> U_new, int N, int t, double g);
double action(std::vector<double> U, int N, double g);

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

  Metropolis_Hasting_Serial_2(U,  moves, N, rd(), a, g, omega, m, epsilon);

  
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

void Metropolis_Hasting_Serial(std::vector<double> &U, int moves, int N, int seed, double a, double g, double omega, double m, double eps){

  std::mt19937 gen(seed);
  std::uniform_real_distribution<double> dis(-a,a);
  std::uniform_real_distribution<double> P(0,1);
  int  t=0; //time slice
  int gap=10;
  double  counter=0.0;
  double X_sq=0.0;

  std::vector <double> U_init; U_init.resize(N);
  std::vector <double> U_new; U_new.resize(N);
  std::fill(U_init.begin(), U_init.end(), 0.0);  
  std::fill(U_new.begin(), U_new.end(), 0.0);
  
  //Variables a calcular
  double ratio=0.0;
  double Delta_S=0.0;  
  double alpha=0;
  
  U_new=U;
  
  for(int jj=0; jj<moves;++jj){
    for(int ii=0; ii<N;++ii){
      t=std::floor(N*P(gen));
      U_init[t]=U[t];
      if(t!=0 && t!=N-1){
	U_new[t]=U_init[t]+dis(gen);
	alpha=P(gen);
	Delta_S=action(U_new,N,g)-action(U,N,g);
	//Delta_S=action_change_2(U_init, U_new, N, t, g);
	ratio=std::min(1.0, std::exp(-Delta_S));
	if(alpha<ratio){
	  U[t]=U_new[t];
	}	
      } else {U[t]=U_init[t];}
    }
    if (jj%gap==0){
      ++counter;
      for (int k=0; k<N; ++k){
	X_sq += std::pow( transform_u_to_x(U[k], m,eps, omega), 2);
      }
    }
  }
  double mean_sq= X_sq/(N*counter);
  std::cout<<"The expected value of x² is "<<mean_sq<<std::endl;
  // return U;
  for(int pp=0;pp<N;++pp){
    std::cout<<pp<<"  "<<U_init[pp]<<"  "<<U[pp]<<"   "<<U_new[pp]<<std::endl;
 }
}

double action_change( std::vector<double> U_init, std::vector<double> U_new, int N, int t, double g){
  
  double ds = 0.0;
  if ( (0!=t) & (N-1 != t) ){
    ds = std::pow(U_new[t], 2) - std::pow(U_init[t], 2) - g*(U_new[t]-U_init[t])*(U_init[t+1] + U_init[t-1]);
  }
  return ds;
}

double transform_u_to_x(double u, const double M, const double EPS, const double OMEGA){
  double k = std::pow(EPS * OMEGA *0.5, 2);
  double x = std::sqrt(EPS/(M*(k+1)))*u;
  return x;
}


double action_change_2( std::vector<double> U_init, std::vector<double> U_new, int N, int t, double g){
  
  double ds = 0.0;
  if ( (0!=t) & (N-1 != t) ){
    ds = std::pow(U_new[t], 2) - std::pow(U_init[t], 2) - g*(U_new[t]*U_init[t+1] - U_init[t]*U_init[t+1]);
  }
  return ds;
}

double action(std::vector<double> U, int N, double g){
  
  double s=0.0;
  double first_term=0.0;
  double second_term=0.0;
  for(int tt=1;tt<N;++tt){
    first_term += U[tt]*U[tt]; 
  }
  for(int jj=0;jj<N-1;++jj){
    second_term += U[jj]*U[jj+1];    
  }
  s=first_term-g*second_term;
  return s;
}

void Metropolis_Hasting_Serial_2(std::vector<double> &U, int Nsweeps, int N, int seed, double a, double g, double omega, double m, double eps){

  std::mt19937 gen(seed);
  std::uniform_real_distribution<double> dis(-a,a);
  std::uniform_real_distribution<double> P(0.0,1.0);
  //int  t=0; //time slice
  int gap=10;
  double  counter=0.0;
  double X_sq=0.0;

  
  //Variables a calcular
  double ratio=0.0;
  double Delta_S=0.0;  
  double alpha=0;
  double shift=0.0;
  
 
  
  for(int jj=1; jj<=Nsweeps;jj++){
    for(int t=1; t<N-1;t++){
      //t=std::floor((N-1)*P(gen));
      //if(t!=0 && t!=N-1){
	shift=dis(gen);
	alpha=P(gen);
	Delta_S=shift*(shift+2.0*U[t]-g*(U[t-1]+U[t+1]));
	//Delta_S=action_change_2(U_init, U_new, N, t, g);
	ratio=std::min(1.0, std::exp(-Delta_S));
	if(alpha<=ratio){
	  U[t] += shift;
	}	
	//} 
    }
    if (jj%gap==0){
      ++counter;
      for (int k=0; k<N; ++k){
	X_sq += std::pow( transform_u_to_x(U[k], m,eps, omega), 2);
      }
      //double mean_sq= X_sq/(N*counter);
	//std::cout<<jj <<"  "<<mean_sq<<std::endl;
    }
  }
   double mean_sq= X_sq/(N*counter);//transform_u_to_x(X_sq,m,eps,omega)/(N*counter);
  std::cout<<"The expected value of x² is "<<mean_sq<<std::endl;
  // return U;
}
