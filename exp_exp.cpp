#include<iostream>
#include<cmath>
#include<vector>
#include<random>

int main(int argc, char **argv){
  //Declaramos las variables del problema
  double mean = 0;
  double mean_sq=0;
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

  //Comienza el sampleo.
  for (int ii = 0; ii<samples; ++ii){

    ++seed;//Variamos la semilla en cada ciclo.

    R = dis(gen); //Establecemos el parámetro aleatorio de la propuesta de actualización.

    x_p = x_val + step*R; // Se propone la catualización.

    DX= x_p*x_p - x_val*x_val; //Se calcula el cambio en el argumento de la exponencial.
    ratio = std::exp(-DX); //Se establece el ratio de aceptación.
    u = P(gen); //Se genera un número aleatorio entre 0 y 1.

    if (u<ratio){ 
      x_val=x_p;  //Actualizamos según ratio de aceptación.
    }


    sum += x_val;
    sum_sq += x_val*x_val; // Actualizamos las variables importantes para calcilar el promedio y el promedio de los cuadrados.
    
    
     }

 
  mean = sum /samples; // Se calcula el promedio.

  sigma = std::sqrt((sum_sq/samples)-mean*mean); // Se calcula la desviación standard.

  mean_sq = sum_sq/samples; // Promedio de los cuadrados.

  std::cout<<mean<<" "<<mean_sq <<" "<< sigma<<std::endl; //Se imprimen los cálculos.


  return 0;







}
