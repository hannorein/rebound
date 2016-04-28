#!/bin/gnuplot
set key top left
set xlabel "time [years]"
set ylabel "semimajor axis [AU]"
set multiplot layout 2,1
beta = 0.01
set lmargin 12
k = 2.497557889905430e-03*beta*4./3.
a(t) = 0.001*sqrt(1.-k*t)
set st d  l
set autoscale xfix 
set xtics 1000

plot "a.txt" u ($1/2./pi):($2) notit

set ylabel "semimajor axis error [AU]"
plot "a.txt" u ($1/2./pi):(abs($2-a($1/2./pi))) notit

pause -1
