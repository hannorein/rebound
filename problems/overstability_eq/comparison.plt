#!/bin/gnuplot

#!./makeheights

set output "comparison.pdf"
set term pdf monochrome dashed enhanced size 8in,4in
set multiplot
set xlabel "optical depth ~@{&'{/Symbol t}}{.4{/Symbol -}}"
set key top left
set st d lp
set autoscale fix
set size 0.33,0.5

set origin 0,0.5
set lmargin at screen 0+0.05
set rmargin at screen 0.33-0.01
set tmargin at screen 0+0.97
set bmargin at screen 0.57
set ylabel "filling factor FF_0"
plot "Kinetic-0.5-1.dat" u 2:1 t "kinetic theory" w l, "<cat heights.txt | sort -g" u 1:($1/$2) t "simulation" ls 1

set origin 0.33,0.5
set lmargin at screen 0.33+0.05
set rmargin at screen 0.66-0.01
set ylabel "viscosity {/Symbol n}"
plot \
'Kinetic-0.5-1.dat' 		u 2:($7) 			t "{/Symbol n}_{col} kinetic theory" w l ls 1, \
'Kinetic-0.5-1.dat' 		u 2:(-2./3.*$5) 		t "{/Symbol n}_{trans} kinetic theory" w l ls 2, \
"<cat out*/data.txt | sort -g" 	u 1:($7) 			t "{/Symbol n}_{col} simulation" ls 3 w p, \
"<cat out*/data.txt | sort -g" 	u 1:(2./3.*$12) 		t "{/Symbol n}_{trans} simulation" ls 4 w p
#"<cat out*/data.txt | sort -g" u 1:($7+2./3.*$12) t "{/Symbol n}_{total} simulation" ls 1, \
#"Kinetic-0.5-1.dat" u 2:($7-2./3.*$5) t "{/Symbol n}_{total} kinetic theory" w l, \

set origin 0.66,0.5
set lmargin at screen 0.66+0.05
set rmargin at screen 0.99
set ylabel "velocity dispersion c^2"
plot "Kinetic-0.5-1.dat" u 2:($3*$3) t "kinetic theory" w l, "<cat out*/data.txt | sort -g" u 1:($6/3.) t "simulation" ls 1



set origin 0,0
set tmargin at screen -0.5+0.97
set bmargin at screen 0.07
set lmargin at screen 0+0.05
set rmargin at screen 0.33-0.01
set ylabel "{/Symbol P}_{xx}"
plot "Kinetic-0.5-1.dat" u 2:4 t "kinetic theory" w l, "<cat out*/data.txt | sort -g" u 1:(1./3.*($10+$11)-2./3.*$9) t "simulation" ls 1

set origin 0.33,0
set lmargin at screen 0.33+0.05
set rmargin at screen 0.66-0.01
set ylabel "{/Symbol P}_{yy}"
plot "Kinetic-0.5-1.dat" u 2:6 t "kinetic theory" w l, "<cat out*/data.txt | sort -g" u 1:(1./3.*($9+$11)-2./3.*$10) t "simulation" ls 1

set origin 0.66,0
set lmargin at screen 0.66+0.05
set rmargin at screen 0.99
set ylabel "{/Symbol P}_{xy}"
plot "Kinetic-0.5-1.dat" u 2:5 t "kinetic theory" w l, "<cat out*/data.txt | sort -g" u 1:(-$12) t "simulation" ls 1



