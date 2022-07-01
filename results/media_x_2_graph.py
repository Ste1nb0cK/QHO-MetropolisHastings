import math
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

def f1(x):
  return 0.5

#lista con todos los pasos de Metropolis usados
Paso = [10, 20, 100, 200, 1000, 2000, 10000, 20000, 100000, 200000]

#empezamos por traer los archivos
df500 = pd.read_csv("datos_x.txt", header=0, sep="\t", usecols=[0,1],names=["moves","x_square_value"])

#Para un valor fijo de paso, extraemos el promedio y desviación estandar del valor de x²
x_square_value_avg_list500 = []
x_square_value_std_list500 = []

for p in Paso:
    #aquí cálculamos el promedio para todos los elementos que tienen el moves dado y sacamos el promedio
    x_square_value_avg500 = df500[df500['moves']== p]["x_square_value"].mean()
    x_square_value_avg_list500.append( x_square_value_avg500)

    #aquí hacemos esencialmente lo mismo pero cálculamos la desviación estandar
    x_square_value_std500 = df500[df500['moves']== p]["x_square_value"].std()
    x_square_value_std_list500.append( x_square_value_std500)

#aquí hacemos las graficas
plt.figure(figsize=(13,10))
l = np.arange(10,200000,100)
plt.plot(l, [f1(i) for i in l], label = "Teorico")
plt.errorbar(Paso, x_square_value_avg_list500, yerr = x_square_value_std_list500, marker = 'o', label = 'N=500')
plt.title('Valor medio de $x^2$ en función del paso de Metropolis en escala logarítmica')
plt.xscale("log")
#plt.yscale("log")
plt.xlabel(r'Pasos de metropolis $S$')
plt.ylabel(r'Promedio de x al cuadrado $E(500,S)$')
plt.legend()
plt.grid()
plt.savefig("Xcuadrado.png")
