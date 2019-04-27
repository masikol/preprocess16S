#!/usr/bin/gnuplot -persist -c

# First argument -- file with data
# Second argument -- output file.

set terminal png enhanced
set output ARG2
set autoscale y
set grid
set xlabel "Avg-quality"
set ylabel "Number-of-reads"
set key autotitle columnhead


plot ARG1 using 1:2 with lines lt rgb 'red' lw 2

