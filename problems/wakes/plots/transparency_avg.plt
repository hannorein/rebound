#!/bin/gnuplot
set terminal pdf enhanced color size 6in,5in
set output "transparency_avg.pdf"
#set palette rgbformulae 30,13,10    #rainbow palette
set palette defined ( 0 '#000090',\
                      1 '#000fff',\
                      2 '#0090ff',\
                      3 '#0fffee',\
                      4 '#90ff70',\
                      5 '#ffee00',\
                      6 '#ff7000',\
                      7 '#ee0000',\
                      8 '#7f0000')
set cblabel "Optical Depth"
#set cbrange [0.55:0.7]
set autoscale cbfix
#unset xtics
#unset ytics
#set format r "%.0f {/Symbol \260}"
#set rtics 15
set autoscale fix
r(B)=(pi/2.-B)*180./pi
psi(B)=.12*r(B)+2
set hidden3d
set pm3d map 
set pm3d corners2color mean
splot "transparency_avg.txt" u (r($1)*cos($2)):(r($1)*sin($2)):(psi($1)):3   notit
set polar
set grid polar
unset hidden3d
unset pm3d
