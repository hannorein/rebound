#!/sw/bin/gnuplot
set terminal pdf enhanced monochrome dashed size 3 in, 2.5 in
set xlabel "number of OMP processes per MPI node"
set logscale y 
set yrange [0.2:100]
set xrange [0.8:8.2]
set ylabel "timesteps per second"
set key top left
set xtics 0,1
set x2tics ("64" 1,"32" 2, "16" 4,"8" 8)
set xtics ("1" 1,"2" 2, "4" 4,"8" 8)
set x2label "number of MPI nodes"
set output "scaling_ompmpi.pdf"
plot "<sort -g -k2 scaling10k.txt" u (64./$2):(100./$3) w lp lt 4 t "OMP+MPI, 64 processes in total, 10k particles","<sort -g -k2 scaling50k.txt" u (64./$2):(100./$3) w lp lt 4 t "OMP+MPI, 64 processes in total, 50k particles","<sort -g -k2 scaling200k.txt" u (64./$2):(100./$3) w lp lt 4 t "OMP+MPI, 64 processes in total, 200k particles","<sort -g -k2 scaling500k.txt" u (64./$2):(100./$3) w lp lt 4 t "OMP+MPI, 64 processes in total, 500k particles"
