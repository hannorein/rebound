#!/usr/bin/gnuplot
set term pdf
set output "growthrate.pdf"
set xlabel "g"
set ylabel "growthrate [Omega]"
set multiplot
plot "killtime.txt" u 1:(11.5129/$2) notit, "henriksguess.txt" t "Henrik's guess" w l
set size 0.4,0.4
set origin 0.05,0.55
unset xlabel
unset ylabel
set xrange [0.095:0.1]
plot "killtime.txt" u 1:(11.5129/$2) notit, "henriksguess.txt" notit w l

