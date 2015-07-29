#!/bin/gnuplot
set key top left
set xlabel "time [years]"
set multiplot layout 2,1
beta = 0.01
set lmargin 12
k = 2.497557889905430e-03*beta
a(t) = sqrt(1.-k*t)
set st d  l
set autoscale xfix 
set xtics 10000

set ytics 0.1
set ylabel "semimajor axis [AU]"
plot "radius.txt" u ($1/2./pi):($2) notit

set ytics 1e-4
set ylabel "semimajor axis error [AU]"
plot "radius.txt" u ($1/2./pi):($2-a($1/2./pi)) notit

pause -1
