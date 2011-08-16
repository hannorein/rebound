#!/sw/bin/gnuplot
set terminal pdf enhanced monochrome dashed size 3 in, 2.3 in
set xlabel "number of nodes"
set logscale 
set ylabel "timesteps per second"
set key top left
set output "scaling_strong.pdf"
a = 0.05
set yrange [0.03:50]
set xrange [0.9:100]
f800k(x) = a800k*x
f400k(x) = a400k*x
f200k(x) = a200k*x
f100k(x) = a100k*x
f50k(x)  = a50k *x
f12k(x)  = a12k *x
fit [0:80] f800k(x) "scaling800k.txt" u 2:(100./$3) via a800k
fit [0:50] f400k(x) "scaling400k.txt" u 2:(100./$3) via a400k
fit [0:50] f200k(x) "scaling200k.txt" u 2:(100./$3) via a200k
fit [0:50] f100k(x) "scaling100k.txt" u 2:(100./$3) via a100k
fit [0:30] f50k(x)  "scaling50k.txt"  u 2:(100./$3) via a50k
fit [0:30] f12k(x)  "scaling12.5k.txt"  u 2:(100./$3) via a12k

plot \
"<sort -g -k2 scaling12.5k.txt" u 2:(100./$3) w lp lt 7  t "MPI, 12.5k particles", \
"<sort -g -k2 scaling50k.txt" u 2:(100./$3) w lp lt 2  t "MPI, 50k particles", \
"<sort -g -k2 scaling200k.txt" u 2:(100./$3) w lp lt 4 t "MPI, 200k particles", \
"<sort -g -k2 scaling800k.txt" u 2:(100./$3) w lp lt 6 t "MPI, 800k particles", \
f200k(x) t "linear scaling" lt 1, \
f50k(x) notit lt 1, \
f800k(x) notit lt 1, \
f12k(x) notit lt 1, \
f200k(x) notit lt 1


#"<sort -g -k2 scaling100k.txt" u 2:(100./$3) w lp lt 3 t "MPI, 100k particles", \
#"<sort -g -k2 scaling400k.txt" u 2:(100./$3) w lp lt 5 t "MPI, 400k particles", \
#f400k(x) notit lt 1
#f100k(x) notit lt 1, \
