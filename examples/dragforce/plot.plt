#!/usr/bin/gnuplot
set xlabel "time"
set multiplot layout 2,1
set ylabel "position"
unset key
set autoscale fix
set lmargin 12
plot "r.txt" u 1:2
set ylabel "absolute position error"
set logscale y
set format y "10^{%T}"
set ytics 100
plot "r.txt" u 1:(abs($2-exp(-$1)))

pause -1
