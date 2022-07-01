#!/bin/bash

#Ciclo para poder obtener los datos del programa
for step in 10 20 100 200 1000 2000 10000 20000 100000 200000
  do
    mpirun -np 6 ./foo 500 $step 1 >> datos.txt
  done
#ejecución de los scripst de python
python3 results/media_x_2_graph.py
