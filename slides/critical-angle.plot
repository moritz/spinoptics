set encoding iso_8859_1
set term postscript landscape font 30
set output 'critical-angle.ps'
set style line 1 linewidth 5 linecolor 1
# set style labels 
set xlabel "rel. alpha"
set ylabel "angle / °"
plot [0:0.20] - 180 / pi * asin(x - sqrt(x*x + 1)) title "Critical angle" ls 1

# vim: set ft=gnuplot fileencoding=latin1
