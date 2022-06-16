 #include<iostream>
#include<cmath>
#include<vector>
#include<random>
#include<cstdlib>

double action_change(const std::vector<double> X, const std::vector<double> X_old, const std::vector<double> X_new, int t,
		     double omega, double m, int T);

int main(int argc, char **argv){

  double m = 2;
  double omega = 1;
  int T = std::atoi(argv[1]);
  double N = std::atoi(argv[2]);
  int seed = 1;
  double counter = 0;
  double X_sq=0;
  
  std::mt19937 gen(seed);
  std::uniform_real_distribution<double> dis(-0.5,0.5);
  std::uniform_real_distribution<double> P(0,1);
  double ratio = 0;
  double DS = 0;
  double u = 0;
  int gap = 10;

  std::vector<double> X;
  X.resize(T);
  std::vector<double> X_old;
  X_old.resize(T);
  std::vector<double> X_new;
  X_new.resize(T);
  std::vector<double> corr_data;
  corr_data.resize(T);
  int t = 0;
  int burn_in_st = 100000;

  for(int ii=0; ii<T; ++ii){
    ++seed;
    X[ii] = dis(gen);
    X_old[ii]=0;
    X_new[ii]=0;
    corr_data[ii]=0;
  }

  for (int ii =0; ii<N; ++ii){
    for (int ii = 0; ii<T; ++ii){
      t = int(T*P(gen));
      X_old[t] = X[t];
      X_new[t] = X_old[t] + dis(gen);
      DS = action_change(X, X_old, X_new, t, omega, m, T);
      u = P(gen);
      ratio = std::exp(-DS);
      if (u<ratio){
	X[t] = X_new[t];
      }
      else { X[t]=X_old[t];}
    }

    if (ii%gap==0){

      ++counter;
      for (int ii=0; ii<T; ++ii){
      X_sq += X[ii]*X[ii];
      corr_data[ii] = X[ii]*X[ii-5];

      }}

    
    
  }

  
  double mean_sq= X_sq/(2*T*counter);
  std::cout<<"The expected value of x² is "<<mean_sq<<std::endl;

  std::cout << "The correlator data is"<< std::endl;
  for (int ii = 5; ii<T; ++ii ){
    std::cout<<corr_data[ii]<<" "<<std::endl;
   
  }

  
  
  return 0;
  }

double action_change(const std::vector<double> X, const std::vector<double> X_old, const std::vector<double> X_new, int t, double omega,
		     double m, int T){

  double DS = 0;
  if (t!=T-1){
  DS =0.5*m*(std::pow(X[t+1]-X_new[t],2)-std::pow(X[t+1]-X_old[t],2) + 0.25*omega*omega*(std::pow(X[t+1]+X_new[t],2)
											-std::pow(X[t+1]+X_old[t],2)));}

   else{
  DS =0.5*m*(std::pow(X[0]-X_new[t],2)-std::pow(X[0]-X_old[t],2) + 0.25*omega*omega*(std::pow(X[0]+X_new[t],2)
  -std::pow(X[0]+X_old[t],2)));}
  return DS;
}
