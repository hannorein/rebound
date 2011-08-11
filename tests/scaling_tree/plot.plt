#!/sw/bin/gnuplot
set terminal png
set xlabel "number of particles"
set logscale 
set ylabel "cpu time [s]"
set key top left

set output "plot.png"

f(x) = a*x*log(x)

fit f(x) "scaling.txt" u 1:3 via a 
plot "scaling.txt" u 1:3 w lp t "single processor", f(x) t "N log(N)"
