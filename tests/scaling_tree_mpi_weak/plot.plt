#!/sw/bin/gnuplot
set terminal pdf enhanced monochrome dashed size 3 in, 2.3 in
set xlabel "number of nodes"
set logscale xy 
set ylabel "timesteps per second"
set key top left
set output "scaling_weak.pdf"
set yrange [0.04:1]
set xrange [0.9:100]
f25k(x) = c25k/log(b25k*x+0.0000001)
f50k(x) = c50k/log(b50k*x+0.0000001)
f100k(x) = c100k/log(b100k*x+0.0000001)
fit f25k(x) "scaling25k.txt" u 2:(100./$3) via b25k,c25k
fit f50k(x) "scaling50k.txt" u 2:(100./$3) via b50k,c50k
fit f100k(x) "scaling100k.txt" u 2:(100./$3) via b100k,c100k

plot "<sort -g -k2 scaling25k.txt" u 2:(100./$3) w lp t "MPI, 25k particles per node" lt 2, "<sort -g -k2 scaling50k.txt" u 2:(100./$3) w lp t "MPI, 50k particles per node" lt 3,  "<sort -g -k2 scaling100k.txt" u 2:(100./$3) w lp t "MPI, 100k particles per node" lt 4 , f25k(x) t "1/log(k)" lt 1, f50k(x) notit lt 1 , f100k(x) notit lt 1


