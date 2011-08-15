#!/sw/bin/gnuplot
set terminal pdf enhanced monochrome dashed size 3 in, 2.3 in
set xlabel "number of particles"
set logscale 
set ylabel "timesteps per second"
set key top right
set output "scaling_collisions.pdf"
set xrange [10:50000]
set yrange [10:1e6]
flinear(x) = alinear/x
fquadratic(x) = aquadratic/(x*x)
fnlogn(x) = anlogn/(x*log(x))
fit flinear(x) "scaling_sweep.txt" u 1:(500./$2) via alinear
fit fquadratic(x) "scaling_direct.txt" u 1:(500./$2) via aquadratic
fit fnlogn(x) "scaling_tree.txt" u 1:(500./$2) via anlogn

plot "<sort -g -k1 scaling_sweep.txt" u 1:(500./$2) w lp lt 2  t "Sweeping algorithm", flinear(x) lt 1 t "1/N","<sort -g -k1 scaling_direct.txt" u 1:(500./$2) w lp lt 2  t "Direct collision search", fquadratic(x) lt 1 t "1/N^2","<sort -g -k1 scaling_tree.txt" u 1:(500./$2) w lp lt 2  t "Nearest neighbor search with tree", fnlogn(x) lt 1 t "1/(N log(N))"

