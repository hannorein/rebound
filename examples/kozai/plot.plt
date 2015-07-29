#!/bin/gnuplot
set key top left
set xlabel "time [{/Symbol W}^{-1}]"
set multiplot
set lmargin 14
set rmargin 4
set size 1,0.5

set ylabel "1-eccentricity"
set origin 0,0.5
set logscale y
set yrange [1e-5:1.3]
plot "< awk '{if(NR%2==1) print $0;}' orbits.txt" u 1:(1-$3) notit w l 

set ylabel "inclination [deg]"
set origin 0,0
unset logscale y
set yrange [*:*]
plot "< awk '{if(NR%2==1) print $0;}' orbits.txt" u 1:($4/pi*180.) notit w l 

pause -1
