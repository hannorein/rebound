#!/sw/bin/gnuplot
set terminal png
set xlabel "opening angle"
set logscale y
set ylabel "average relative acceleration error"

set output "plot1.png"
plot "error-mono.txt" u 13:3 w lp t "monopole",  'error-quad.txt' u 13:3 w lp t "quadrupole" 

#set output "plot2.png"
#set xlabel "average relative acceleration error"
#set ylabel "runtime [s]"
#set logscale xy
#plot "error-mono.txt" u 3:6 w lp t "monopole",  'error-quad.txt' u 3:6 w lp t "quadrupole" 
