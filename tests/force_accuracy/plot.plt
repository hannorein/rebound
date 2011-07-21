#!/sw/bin/gnuplot
set output "plot.png"
set terminal png
set xlabel "opening angle"
set ylabel "log(average relative acceleration error) in \%"

plot "error-mono.txt" u 13:(log10($3*100)) w lp t "monopole",  'error-quad.txt' u 13:(log10($3*100)) w lp t "quadrupole" 
