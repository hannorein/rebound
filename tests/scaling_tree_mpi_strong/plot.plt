#!/sw/bin/gnuplot
set terminal png
set xlabel "number of nodes"
set logscale 
set ylabel "timesteps per second"
set key top left
set output "plot.png"
a = 0.05
set yrange [0.05:]
f(x) = a*x
fit [0:50] f(x) "scaling.txt" u 2:(100./$3) via a
plot "scaling.txt" u 2:(100./$3) w lp t "MPI, 1 process per node, 100k particles", f(x) t "linear scaling"

