import math
import matplotlib.pyplot as plt
import pandas as pd

#lista con todos los pasos de Metropolis usados
Paso = [10, 20, 100, 200, 1000, 2000, 10000, 20000, 100000, 200000, 1000000, 2000000, 10000000, 20000000]

#empezamos por traer los archivos
df100 = pd.read_csv("datos100.txt", sep="\t", usecols=[1,2,3,4,5],names=["N","moves","omega","x_square_value", "time"])
df500 = pd.read_csv("datos500.txt", sep="\t", usecols=[1,2,3,4,5],names=["N","moves","omega","x_square_value", "time"])
df1000 = pd.read_csv("datos1000.txt", sep="\t", usecols=[1,2,3,4,5],names=["N","moves","omega","x_square_value","time"])
df1500 = pd.read_csv("datos1500.txt", sep="\t", usecols=[1,2,3,4,5],names=["N","moves","omega","x_square_value", "time"])
df2000= pd.read_csv("datos2000.txt", sep="\t", usecols=[1,2,3,4,5],names=["N","moves","omega","x_square_value", "time"])
#Para un valor fijo de paso, extraemos el promedio y desviación estandar del valor de x²
x_square_value_avg_list100 = []
x_square_value_std_list100 = []

x_square_value_avg_list500 = []
x_square_value_std_list500 = []

x_square_value_avg_list1000 = []
x_square_value_std_list1000 = []

x_square_value_avg_list1500 = []
x_square_value_std_list1500 = []

x_square_value_avg_list2000 = []
x_square_value_std_list2000 = []

for p in Paso:
    #aquí cálculamos el promedio para todos los elementos que tienen el moves dado y sacamos el promedio
    x_square_value_avg100 = df100[df100['p']== p]["x_square_value"].mean()
    x_square_value_avg_list100.append( x_square_value_avg100)

    x_square_value_avg500 = df500[df500['p']== p]["x_square_value"].mean()
    x_square_value_avg_list500.append( x_square_value_avg500)

    x_square_value_avg1000 = df1000[df1000['p']== p]["x_square_value"].mean()
    x_square_value_avg_list1000.append( x_square_value_avg1000)

    x_square_value_avg1500 = df1500[df1500['p']== p]["x_square_value"].mean()
    x_square_value_avg_list1500.append( x_square_value_avg1500)

    x_square_value_avg2000 = df2000[df2000['p']== p]["x_square_value"].mean()
    x_square_value_avg_list2000.append( x_square_value_avg2000)

    #aquí hacemos esencialmente lo mismo pero cálculamos la desviación estandar
    x_square_value_std100 = df100[df100['p']== p]["x_square_value"].std()
    x_sqaure_value_std_list100.append( x_square_value_std100)

    x_sqaure_value_std500 = df500[df500['p']== p]["x_square_value"].std()
    x_square_value_std_list500.append( x_square_value_std500)

    x_square_value_std1000 = df1000[df1000['p']== p]["x_square_value"].std()
    x_square_value_std_list1000.append( x_square_value_std1000)

    x_square_value_std1500 = df1500[df1500['p']== p]["x_square_value"].std()
    x_square_value_std_list1500.append( x_square_value_std1500)

    x_square_value_std2000 = df2000[df2000['p']== p]["x_square_value"].std()
    x_square_value_std_list2000.append( x_square_value_std2000)

#aquí hacemos las graficas
plt.figure(figsize=(13,10))
plt.errorbar(Paso, x_square_value_avg_list100, yerr = x_square_value_std_list100, marker = 'o', label = 'N=100')
plt.errorbar(Paso, x_square_value_avg_list500, yerr = x_square_value_std_list500, marker = 'o', label = 'N=500')
plt.errorbar(Paso, x_square_value_avg_list1000, yerr = x_square_value_std_list1000, marker = 'o', label = 'N=1000')
plt.errorbar(Paso, x_square_value_avg_list1500, yerr = x_square_value_std_list1500, marker = 'o', label = 'N=1500')
plt.errorbar(Paso, x_square_value_avg_list2000, yerr = x_square_value_std_list2000, marker = 'o', label = 'N=2000')
plt.title('Valor medio de $x^2$ en función del paso de Metropolis')
plt.xlabel(r'Pasos de metropolis $S$')
plt.ylabel(r'Promedio de x al cuadrado $E(N,S)$')
plt.legend()
plt.grid()
plt.savefig("Xcuadrado.png")
