#!/bin/bash

#Ciclo para poder obtener los datos del programa variando el paso de Metropolis
for step in 10 20 100 200 1000 2000 10000 20000 100000 200000
do
  mpirun -np 6 ./foo 500 $step 1 >> datos_x.txt
done

#Ejecución de los scripst de python para el paso
python3 results/media_x_2_graph.py

#Ejecucion de los scripts de gnuplot
gnuplot results/speed.p
gnuplot results/parallel.p

#Ciclo para poder obtener llos datos del programa variando m
for m in 1 2 3 4 5 6 7 8
do
  mpirun -np 6 ./foo 500 10000 $m >> datos_m.txt
done

#Ejecución del script de python para la masa
python3 results/oscilador.py
