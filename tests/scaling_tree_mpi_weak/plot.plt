#!/sw/bin/gnuplot
set terminal png
set xlabel "number of nodes"
set logscale x 
set ylabel "timesteps per second"
set key top left
set output "plot.png"
a = 0.05
set yrange [0.:0.2]
f(x) = a
fit f(x) "scaling.txt" u 2:(100./$3) via a
plot "scaling.txt" u 2:(100./$3) w lp t "MPI, 50k particles per node", f(x) t "constant scaling"

