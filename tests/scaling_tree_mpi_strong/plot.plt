#!/sw/bin/gnuplot
set terminal png
set xlabel "number of cores"
set logscale 
set ylabel "cpu time [s]"
set key top left

set output "plot.png"


plot "scaling.txt" u 2:3 w lp 
