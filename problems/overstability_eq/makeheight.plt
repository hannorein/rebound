#!/bin/gnuplot
N=system("wc -l heights.txt | awk '{print $1}'")
H=5


f(z) = (z<-H?0:N*((z>H?1.:(2./3.*H+z-z**3/3./H**2)/H/4.*3.)))
fit f(x) "<pgt_comm heights.txt 0" via H

set table "height.txt"
plot H
unset table
