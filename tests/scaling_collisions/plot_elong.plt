#!/sw/bin/gnuplot
set terminal pdf enhanced monochrome dashed size 3 in, 2.3 in
set xlabel "number of particles N"
set logscale 
set ylabel "timesteps per second"
set key top right
set output "scaling_collisions_elong.pdf"
set xrange [10:20000]
set yrange [10:100000]
flinear(x) = alinear/x
fquadratic(x) = aquadratic/(x*x)
fnlogn(x) = anlogn/(x*log(x))
fit [100:] flinear(x) "scaling_sweep_elong.txt" u 1:(500./$2) via alinear
fit [100:] fquadratic(x) "scaling_direct_elong.txt" u 1:(500./$2) via aquadratic
fit [100:]fnlogn(x) "scaling_tree_elong.txt" u 1:(500./$2) via anlogn

plot "<sort -g -k1 scaling_sweep_elong.txt" u 1:(500./$2) w lp lt 2  t "plane-sweep algorithm", flinear(x) lt 1 t "1/N", "<sort -g -k1 scaling_direct_elong.txt" u 1:(500./$2) w lp lt 3  t "direct collision search", "<sort -g -k1 scaling_tree_elong.txt" u 1:(500./$2) w lp lt 4  t "nearest neighbor search with tree", fquadratic(x) lt 1 t "1/N^2"

