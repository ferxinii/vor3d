#!/bin/bash

gnuplot <<- EOF
        set terminal pngcairo enhanced font 'Arial,18' size 1080,1080 enhanced

        set output "$1_0.png"

        binwidth = 0.5;
        bin(x,width)=width*floor(x/width) + width/2.0

        set title "Distribution of cell volumes"
        set xlabel 'Volume (cm^3)'
        set ylabel 'Count'

        set xtics 2
        set grid
        set boxwidth binwidth
        set style fill solid border -1 
        set xrange[0:60]
        plot '$1.txt' using (bin(\$5,binwidth)):(1.0) smooth freq with boxes lc rgb "#0072bd" notitle



        set output "$1_1.png"

        set xrange[-26:-2]
        set yrange[0:10]

        plot '$1.txt' using 4:5 with points notitle
EOF









