#!/bin/bash

#Ciclo para poder obtener los datos del programa
for N in 100 500 1000 1500 2000
do
  for step in 10 20 100 200 1000 2000 10000 20000 100000 200000 1000000 2000000 10000000 20000000
  do
    ./a.out $N $step 1 >> datos$N.txt
  done
done

#ejecución de los scripst de python
python3 results/media_x_2_graph.py
