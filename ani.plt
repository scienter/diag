#t=200000
unset multiplot

#set terminal pngcairo size 500,400 enhanced font 'Verdana,10' enhanced
set terminal pngcairo size 800,500 enhanced font 'Verdana,10'

refile=sprintf("summary%06d.png",t)
set output refile

set multiplot
set view map

dx=4.0e-2
maxE=600
end_time=200000
n0=7e24
rangeX=80
rangeR=50
a0_max=12

unset y2tics
unset y2label
set xrange[t*dx-rangeX:t*dx]


set lmargin at screen 0.11
set rmargin at screen 0.40
set bmargin at screen 0.90
set tmargin at screen 0.55

unset xlabel
set yrange [-rangeR:rangeR]
set cbrange [0:10]
set cbtics 2,2,10 font "Times-Roman,15"
unset xtics
set ytics -rangeR,rangeR/2,rangeR nomirror font "Times-Roman,15"
set ylabel 'x [um]' font "Times-Roman,15" offset -2,0
set cblabel 'n0 [norm.]' font "Times-Roman,15" offset 1,0
file=sprintf("0density%d_sum_0",t)
set palette defined (0 1 1 1,0.25 0 0 1,0.5 0 1 0,0.75 1 0 0,1 1 1 0)
spl file u ($1*1e6):($2*1e6):($3/n0) title "z-x plane" w pm


set bmargin at screen 0.2
set tmargin at screen 0.55

unset xtics
set yrange [-rangeR:rangeR]
set cbrange [-5:5]
set cbtics -5,5,5 font "Times-Roman,15"
set cblabel 'a0 [norm.]' font "Times-Roman,15" offset 1,0
set xtics 0,30,8000 nomirror font "Times-Roman,15"
set ytics -rangeR,rangeR/2,rangeR nomirror font "Times-Roman,15"
set ylabel 'x [um]' font "Times-Roman,15" offset -2,0
set xlabel 'z [um]' font "Times-Roman,15" offset 0,-0
file=sprintf("fieldSplit%d_sum_0",t)
set palette defined (0 0 0 1,0.5 1 1 1,1 1 0 0)
spl file u ($1*1e6):($2*1e6):($4) notitle w pm


set lmargin at screen 0.63
set rmargin at screen 0.88
set bmargin at screen 0.62
set tmargin at screen 0.90

unset key
set xlabel 'z [um]' font "Times-Roman,15" offset 0,-0
set ylabel 'energy [MeV]' font "Times-Roman,15" offset -3,0
set cblabel 'charge [a.u.]' font "Times-Roman,15" offset 2,0
set ytics 0,200,maxE nomirror font "Times-Roman,15"
set cbrange [0:10]
set cbtics 0,10,10 font "Times-Roman,15"
set palette defined (0 1 1 1,0.25 0 0 1,0.5 0 1 0,0.75 1 0 0,1 1 1 0)
file=sprintf("denSpec%4d",t)
set yrange [0:maxE]
spl file u ($1*1e6):($2):($3*1.602e-7*1e-6) w pm


unset xtics
unset ytics
unset xlabel
unset ylabel
unset border
set xrange[t*dx-rangeX:t*dx]
set y2range [-0.2:0.4]
unset y2tics
file=sprintf("cen%d_0",t)
pl file u ($1*1e6):2 axes x1y2 w l lt 5 lw 0.1 lc rgb "black"
unset y2range

set lmargin at screen 0.63
set rmargin at screen 0.69

set xrange [0:10]
set yrange [0:maxE]
file=sprintf("spectrum%d",t)
pl file u ($2*1.602e-7):1 w l lw 1 lc rgb "red"
unset arrow





set lmargin at screen 0.63
set rmargin at screen 0.9
set bmargin at screen 0.2
set tmargin at screen 0.47

set xrange [0:8]
set yrange [0:n0*1.1]
set border
set xlabel 'z [mm]' font "Times-Roman,15" offset 0,0
set xtics -2,2,10 nomirror font "Times-Roman,15"
set ylabel 'n0 [SI]' font "Times-Roman,15" offset 5,-0.5
set ytics 0,1.0*n0,1.0*n0 nomirror font "Times-Roman,15"
set arrow 1 from t*dx*1e-3,0 to t*dx*1e-3,1.1*n0 nohead lt 7 lc rgb "blue"
pl 'densityProfile' u ($1*1e3):($2) w l lw 1 lt 2 lc rgb "black"
unset arrow
set xrange [0:16]
set y2range [0:a0_max]
set y2tics 0,2,20 font "Times-Roman,15"
set y2label 'max a0' font "Times-Roman,15" offset 1,0
unset xtics
unset xlabel
unset ylabel
unset ytics
pl 'laserMax' u ($1*dx*1e-3):($2) axes x1y2 w l lw 1 lt 2 lc rgb "red"
unset y2range
unset y2tics

unset multiplot
unset output

set term x11

print "step",t," is done."
t = t + 2000

if(t<end_time) reread;

