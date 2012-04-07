#!/bin/gnuplot
set output "plot.pdf"
set key top left
set ylabel "1-eccentricity"
set xlabel "time [{/Symbol W}^{-1}]"
set terminal pdf monochrome dashed enhanced
set logscale y
set yrange [1e-6:1]
plot "orbits.txt" u 1:3 t "simulation"
