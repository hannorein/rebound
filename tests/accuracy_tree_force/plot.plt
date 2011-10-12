#!/sw/bin/gnuplot
set terminal pdf monochrome dashed size 3in, 2.3in 
set xlabel "opening angle"
set logscale y
set ylabel "average relative acceleration error"
set key top left
set st d lp
set output "accuracy_force.pdf"

plot "error-mono.txt" u 13:3  t "monopole",  'error-quad.txt' u 13:3  t "quadrupole" 
