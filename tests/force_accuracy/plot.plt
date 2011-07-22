#!/sw/bin/gnuplot
set output "plot.png"
set terminal png
set xlabel "opening angle"
set logscale y
set ylabel "average relative acceleration error"

plot "error-mono.txt" u 13:3 w lp t "monopole",  'error-quad.txt' u 13:3 w lp t "quadrupole" 
