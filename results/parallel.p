set terminal png size 500,500

set output "Parallel-efficiency.png"

set title "P vs n / Computador con 6 threads"

set xlabel "Número de threads (n)"

set ylabel "Paralle-efficiency (P)"

p [:7][:] 1 lt -1 t "Teorica", "datos-p.dat" u 1:2 w lp t "Computacional"
