#define CATCH_CONFIG_MAIN
#include "catch2/catch.hpp"
#include "x_sq_Metropolis_Hasting_Serial.hpp"
#include "action_change.hpp"
#include "transform_u_to_x.hpp"
#include <cmath>


TEST_CASE("Mean of X² computed", "[x_sq_Metropolis_Hasting_Serial]"){

  const int  N= 300; //Discretización del tiempo
  const int moves= 10000; //Cantidad de movimientos en el sampleo
  const double epsilon=0.25;
  const double omega=1;
  const double k=std::pow(epsilon*omega*0.5,2);
  const double g=(1-k)/(1+k);
  const double m=1;
  const double a = 0.5;
  const double seed=1;
  std::vector<double> U;
  U.resize(N);
  
  //Se testea el el algoritmo para calcular <x²> para diferentes masas, y se coompara con el valor teórico
  double x_sq1= x_sq_Metropolis_Hasting_Serial(U, moves, N, seed, a, g, omega, 1, epsilon);
  double x_sq2= x_sq_Metropolis_Hasting_Serial(U, moves, N, seed, a, g, omega, 2, epsilon);
  double x_sq3= x_sq_Metropolis_Hasting_Serial(U, moves, N, seed, a, g, omega, 3, epsilon);
  double x_sq4= x_sq_Metropolis_Hasting_Serial(U, moves, N, seed, a, g, omega, 4, epsilon);
  

  REQUIRE(std::fabs(x_sq1-0.5)<0.05);
  REQUIRE(std::fabs(x_sq2-0.25)<0.05);
  REQUIRE(std::fabs(x_sq3-0.16)<0.05);
  REQUIRE(std::fabs(x_sq4-0.125)<0.05);



}
