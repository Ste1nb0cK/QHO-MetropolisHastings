#include "Metropolis_Hasting_Parallel.hpp"
#include "action_change.hpp"
#include "transform_u_to_x.hpp"

void Metropolis_Hasting_Parallel(std::vector<double> &U, int Nsweeps, int N, int seed, double a, double g, double omega, double m, double eps, int pid, int np){

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

  int Nt=N/np;

  //for all process
  for(int jj=1; jj<=Nsweeps;jj++){
    for(t=1; t<Nt-1;t++){         //medición tomando cada uno de los t's
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
    if (jj%gap==0){
      ++counter;
      for (int k=0; k<Nt; ++k){
	X_sq += std::pow( transform_u_to_x(U[k], m,eps, omega), 2);
      }
      //double mean_sq= X_sq/(N*counter);
	//std::cout<<jj <<"  "<<mean_sq<<std::endl;
    }
  }
  double mean_sq= X_sq/(Nt*counter);
  double mean_sq_sum;
  MPI_Reduce(&mean_sq, &mean_sq_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  if(0==pid){
    double x_sq=mean_sq_sum/np;
    std::cout<<"The expected value of x² is: "<<x_sq<<std::endl;
  }
}
