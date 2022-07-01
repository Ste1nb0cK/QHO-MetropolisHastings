import math
from matplotlib import pyplot as plt
import pandas as pd
import numpy as np

def f1(x):
  return 1/(2*x)

#lista con todos los pasos de Metropolis usados
M = [1, 2, 3, 4, 5, 6, 7, 8]

#empezamos por traer los archivos
df1 = pd.read_csv("datos_m.txt", header=0, sep="\t", usecols=[1,2],names=["x_square_value","mass"])

#Para un valor fijo de masa, extraemos el promedio y desviación estandar del valor de x²
x_square_value_avg_list1 = []
x_square_value_std_list1 = []

for m in M:
    #aquí cálculamos el promedio para todos los elementos que tienen el masa dada y sacamos el promedio
    x_square_value_avg1 = df1[df1["mass"]== m]["x_square_value"].mean()
    x_square_value_avg_list1.append( x_square_value_avg1)

    #aquí hacemos esencialmente lo mismo pero cálculamos la desviación estandar
    x_square_value_std1 = df1[df1["mass"]== m]["x_square_value"].std()
    x_square_value_std_list1.append( x_square_value_std1)

#aquí hacemos las graficas
plt.figure(figsize=(13,10))
l = np.arange(1, 8 , 0.05)
plt.plot(l, [f1(i) for i in l], label = 'Teoríca')
plt.errorbar(M, x_square_value_avg_list1, yerr = x_square_value_std_list1, marker = 'o', label = 'N=500,Paso=10000')
plt.title('Valor medio de $x^2$ en función de la masa')
#plt.xscale("log")
#plt.yscale("log")
plt.xlabel(r'Masa $m$')
plt.ylabel(r'Promedio de x al cuadrado $E(500,10000,m)$')
plt.legend()
plt.grid()
plt.savefig("Xcuadrad_m.png")
