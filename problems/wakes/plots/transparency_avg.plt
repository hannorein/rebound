#!/bin/gnuplot
set terminal pdf enhanced color size 6in,5in
set output "transparency_avg.pdf"
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
set autoscale cbfix
unset xtics
unset ytics
set xrange [-75:75]
set yrange [-75:75]
unset border 
set autoscale fix
r(B)=(pi/2.-B)*180./pi
psi(B)=.12*r(B)+2
set hidden3d
set pm3d map 
set pm3d corners2color mean


len = 70.

set arrow from len*sin(0.*pi/5),len*cos(0.*pi/5) to -len*sin(0.*pi/5),-len*cos(0.*pi/5) nohead lc rgb 'gray' front
set arrow from len*sin(1.*pi/5),len*cos(1.*pi/5) to -len*sin(1.*pi/5),-len*cos(1.*pi/5) nohead lc rgb 'gray' front
set arrow from len*sin(2.*pi/5),len*cos(2.*pi/5) to -len*sin(2.*pi/5),-len*cos(2.*pi/5) nohead lc rgb 'gray' front
set arrow from len*sin(3.*pi/5),len*cos(3.*pi/5) to -len*sin(3.*pi/5),-len*cos(3.*pi/5) nohead lc rgb 'gray' front
set arrow from len*sin(4.*pi/5),len*cos(4.*pi/5) to -len*sin(4.*pi/5),-len*cos(4.*pi/5) nohead lc rgb 'gray' front

len = 71.
set label "{/Symbol f} = 0{/Symbol \260}"   at len*sin(0.*pi/5),len*cos(0.*pi/5) front
set label "{/Symbol f} = 36{/Symbol \260}"  at len*sin(1.*pi/5),len*cos(1.*pi/5) front
set label "{/Symbol f} = 72{/Symbol \260}"  at len*sin(2.*pi/5),len*cos(2.*pi/5) front
set label "{/Symbol f} = 108{/Symbol \260}" at len*sin(3.*pi/5),len*cos(3.*pi/5) front
set label "{/Symbol f} = 144{/Symbol \260}" at len*sin(4.*pi/5),len*cos(4.*pi/5) front
set label "{/Symbol f} = 180{/Symbol \260}" at len*sin(5.*pi/5),len*cos(5.*pi/5) front
len = 78.
set label "{/Symbol f} = 216{/Symbol \260}" at len*sin(6.*pi/5),len*cos(6.*pi/5) front
set label "{/Symbol f} = 252{/Symbol \260}" at len*sin(7.*pi/5),len*cos(7.*pi/5) front
set label "{/Symbol f} = 288{/Symbol \260}" at len*sin(8.*pi/5),len*cos(8.*pi/5) front
set label "{/Symbol f} = 324{/Symbol \260}" at len*sin(9.*pi/5),len*cos(9.*pi/5) front

splot "transparency_avg.txt" u (r($1)*cos($2)):(r($1)*sin($2)):(psi($1)):3   notit
