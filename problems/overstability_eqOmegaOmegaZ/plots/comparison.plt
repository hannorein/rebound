#!/bin/gnuplot
set output "comparison.pdf"
set term pdf monochrome enhanced size 8in,4in
set multiplot
set xlabel "optical depth {/Symbol t}"
set key top left
set st d lp
set autoscale fix

set size 0.33,0.5

set origin 0,0.5
set ylabel "filling factor FF_0"
plot "Kinetic-0.5-1.dat" u 2:1 t "kinetic theory" w l, "<cat ../heights.txt | sort -g" u 1:($1/$2) t "simulation"

set origin 0.33,0.5
set ylabel "collision frequency"
plot "Kinetic-0.5-1.dat" u 2:8 t "kinetic theory" w l, "<cat ../out*/data.txt | sort -g" u 1:($8) t "simulation"

set origin 0.66,0.5
set ylabel "C^2"
plot "Kinetic-0.5-1.dat" u 2:3 t "kinetic theory" w l, "<cat ../out*/data.txt | sort -g" u 1:($6/3.) t "simulation"


set origin 0,0
set ylabel "{/Symbol P}_{xx}"
plot "Kinetic-0.5-1.dat" u 2:4 t "kinetic theory" w l, "<cat ../out*/data.txt | sort -g" u 1:(1./3.*($10+$11)-2./3.*$9) t "simulation"

set origin 0.33,0
set ylabel "{/Symbol P}_{yy}"
plot "Kinetic-0.5-1.dat" u 2:6 t "kinetic theory" w l, "<cat ../out*/data.txt | sort -g" u 1:(1./3.*($9+$11)-2./3.*$10) t "simulation"

set origin 0.66,0
set ylabel "{/Symbol P}_{xy}"
plot "Kinetic-0.5-1.dat" u 2:5 t "kinetic theory" w l, "<cat ../out*/data.txt | sort -g" u 1:(-$12) t "simulation"

