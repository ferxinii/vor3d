#!/bin/bash

gnuplot <<- EOF
        # ---------------------- VOLUME DISTRIBUTION ------------------------ 
        set terminal pngcairo enhanced font 'Arial,18' size 1080,1080 enhanced

        # --- STATS ---
        unset xrange
        unset yrange
        stats "$1.txt" u 5 name "V" nooutput

        # --- BIN FUN ---
        # binwidth = 0.2;
        binwidth = (ceil(V_max) - floor(V_min))/30.0;
        bin(x,width)=width*floor(x/width) + width/2.0


        set output "$1_0.png"

        set title "Distribution of cell volumes"
        set xlabel 'Volume (cm^3)'
        set ylabel 'Frequency'

        # set xtics 1
        set grid
        set boxwidth binwidth
        set style fill solid border -1 
        set xrange[0:ceil(V_max)]
        plot '$1.txt' using (bin(\$5,binwidth)):(1.0/V_records) smooth freq with boxes lc rgb "#0072bd" notitle



        # ---------------------- VOLUME VS Z ------------------------ 

         # --- STATS ---
        unset xrange
        unset yrange
        stats '$1.txt' u 4 nooutput name "Z"

        # --- FIT ---
        f(x) = m*x + b
        m = 1      # initial slope
        b = 1      # initial intercept
        set fit quiet
        set fit logfile '/dev/null'
        fit f(x) "$1.txt" using 4:5 via m,b 

        set output "$1_1.png"

        set xtics 2
        # set xrange[-26:-2]
        # set yrange[0:8]
        set title "Volume and height"
        set xlabel 'z coordinate of seed (cm)'
        set ylabel 'Volume of Voronoi cell (cm^3)'

        plot '$1.txt' using 4:5 with points lt 6 notitle, \
        [Z_min:Z_max] f(x) w l lt 1 lw 4 lc rgb "red" title "fit"
EOF









