#!/sw/bin/gnuplot
set terminal pdf enhanced monochrome dashed size 3 in, 2.3 in
set xlabel "number of particles N"
set logscale 
set ylabel "timesteps per second"
set key top right
set output "scaling_collisions.pdf"
set xrange [100:1e7]
set yrange [0.01:1e5]
flinear(x) = alinear/x**(1.5)
fquadratic(x) = aquadratic/(x*x)
fnlogn(x) = anlogn/(x*log(x))
fit [2000:] flinear(x) "scaling_sweep.txt" u 1:(10./$2) via alinear
fit fquadratic(x) "scaling_direct.txt" u 1:(10./$2) via aquadratic
fit [2000:]fnlogn(x) "scaling_tree.txt" u 1:(10./$2) via anlogn

plot "<sort -g -k1 scaling_sweep.txt" u 1:(10./$2) w lp lt 2  t "plane-sweep algorithm",flinear(x) lt 1 t "1/N^{1.5}", "<sort -g -k1 scaling_direct.txt" u 1:(10./$2) w lp lt 3  t "direct collision search",fquadratic(x) lt 1 t "1/N^{2}", "<sort -g -k1 scaling_tree.txt" u 1:(10./$2) w lp lt 4  t "nearest neighbor search with tree",  fnlogn(x) lt 1 t "1/(N log(N))"

