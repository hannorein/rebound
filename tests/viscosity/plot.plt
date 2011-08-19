#!/sw/bin/gnuplot
set terminal pdf enhanced monochrome dashed size 3in,3in
set output "viscosity.pdf"
set xrange [0.45:1.1]
set yrange [1e-3:1e2]
set logscale xy
set key top left
set st d lp
plot "viscosity.txt" u 2:3 ls 6 lt 1 t "trans REBOUND", '' u 2:4 ls 8 lt 1 t "coll REBOUND", '' u 2:5 ls 4 lt 1 t "grav REBOUND", "daisaka.txt" u 1:2 ls 6 lt 2 t "trans Daisaka", '' u 1:3 ls 8 lt 2 t "coll Daisaka", '' u 1:4 ls 4 lt 2 t "grav Daisaka"
set xlabel "r_h^*"
set ylabel "{/Symbol n}"
