#include<iostream>
#include<cmath>
#include<vector>
#include<random>

int main(int argc, char **argv){

  double mean = 0;
  double sigma = 0.0;
  double x_val = 1;
  double x_p = 0.0;
  double x_sq = 0.0;
  double sum_sq = 0.0;
  double DX = 0.0;
  double ratio =0;
  double u = 0;
  double samples = std::atof(argv[1]);
  double step = std::atof(argv[2]);
  int seed = std::atoi(argv[3]);
  double a = std::atof(argv[4]);
  double R = 0;
  double sum = 0;
  std::mt19937 gen(seed);
  std::uniform_real_distribution<double> dis(-a,a);
  std::uniform_real_distribution<double> P(0,1);

  for (int ii = 0; ii<samples; ++ii){

    ++seed;
    R = dis(gen);
    x_p = x_val + step*R;

    DX= x_p*x_p - x_val*x_val;
    ratio = std::exp(-DX);
    u = P(gen);

    if (u<ratio){
      x_val=x_p;
    }


    sum += x_val;
    sum_sq += x_val*x_val;
    
    
     }

 
  mean = sum /samples;
  sigma = std::sqrt((sum_sq/samples)-mean*mean);
  double mean_sq = sum_sq/samples;

  std::cout<<mean<<" "<<mean_sq<<std::endl;


  return 0;







}
