set terminal png size 500,500

set output "Speedup.png"

set title "S vs n / Computador con 6 threads"

set xlabel "Número de threads (n)"

set ylabel "Speedup (S)"

p [:7][:] x lt -1 t "Teorica", "datos-s.dat" u 1:2 w lp t "Computacional"
