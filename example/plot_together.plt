#!/bin/bash


# ------- MERGE ALL TXT FILES INTO A SINGLE ONE ------
OUT="together.txt"

# : > "$OUT"                              # start with an empty output
# for f in $(ls -1t -- *.txt); do                      # loop through all .txt in the cwd
#     [[ "$f" == "$OUT" ]] && continue    # skip the output itself if it lives here
#     cat "$f" >> "$OUT"                  # add exactly two blank lines after each file
#     echo >> "$OUT"
#     echo >> "$OUT"
# done

: > "$OUT"                              # start with an empty output
cat "L_adult.txt" >> "$OUT"
cat "R_adult.txt" >> "$OUT"
echo >> "$OUT"
echo >> "$OUT"
cat "L_13yo.txt" >> "$OUT"
cat "R_13yo.txt" >> "$OUT"
echo >> "$OUT"
echo >> "$OUT"
cat "L_8yo.txt" >> "$OUT"
cat "R_8yo.txt" >> "$OUT"
echo >> "$OUT"
echo >> "$OUT"
cat "L_3yo.txt" >> "$OUT"
cat "R_3yo.txt" >> "$OUT"
echo >> "$OUT"
echo >> "$OUT"


gnuplot <<- EOF
        # ---------------------- VOLUME DISTRIBUTION ------------------------ 
        set terminal pdfcairo enhanced font 'Arial,12' size 4,4 enhanced transparent

        set lmargin at screen 0.15
        set rmargin at screen 0.85
        set bmargin at screen 0.15
        set tmargin at screen 0.85
        set size 1,1   # force square canvas

        # --- STATS ---
        unset xrange
        unset yrange
        stats "$OUT" u 5 name "V" nooutput

        stats "$OUT" u 5 i 0 name "V0" nooutput
        stats "$OUT" u 5 i 1 name "V1" nooutput
        stats "$OUT" u 5 i 2 name "V2" nooutput
        stats "$OUT" u 5 i 3 name "V3" nooutput
        print V0_mean
        print V1_mean
        print V2_mean
        print V3_mean
        print V0_records
        print V1_records
        print V2_records
        print V3_records

        # --- BIN FUN ---
        # binwidth = 0.2;
        binwidth = (ceil(V_max) - floor(V_min))/50.0;
        bin(x,width)=width*floor(x/width) + width/2.0


        set output "together_0.pdf"

        set title "Distribution of cell volumes"
        set xlabel 'Volume (cm^3)'
        set ylabel 'Frequency'

        # set xtics 1
        set grid
        set boxwidth binwidth
        set style fill empty border -1
        set xrange[0:ceil(V_max)]
        plot '$OUT' using (bin(\$5,binwidth)):(1.0/V_records) i 0 smooth freq with boxes lw 3 lc rgb "#0072BD" title 'Adult' , \
             '$OUT' using (bin(\$5,binwidth)):(1.0/V_records) i 1 smooth freq with boxes lw 3 lc rgb "#D95319" title '13 year old' , \
             '$OUT' using (bin(\$5,binwidth)):(1.0/V_records) i 2 smooth freq with boxes lw 3 lc rgb "#EDB120" title '8 year old', \
             '$OUT' using (bin(\$5,binwidth)):(1.0/V_records) i 3 smooth freq with boxes lw 3 lc rgb "#7E2F8E" title '3 year old'



        # ---------------------- VOLUME VS Z ------------------------ 

         # --- STATS ---
        unset xrange
        unset yrange
        stats '$OUT' u 4 i 0 nooutput name "Z0"
        stats '$OUT' u 4 i 1 nooutput name "Z1"
        stats '$OUT' u 4 i 2 nooutput name "Z2"
        stats '$OUT' u 4 i 3 nooutput name "Z3"

        # --- FIT ---
        f0(x) = m0*x + b0; m0 = 1; b0 = 1
        f1(x) = m1*x + b1; m1 = 1; b1 = 1
        f2(x) = m2*x + b2; m2 = 1; b2 = 1
        f3(x) = m3*x + b3; m3 = 1; b3 = 1

        set fit quiet
        set fit logfile '/dev/null'

        fit f0(x) "$1" using 4:5 i 0 via m0,b0 
        fit f1(x) "$1" using 4:5 i 1 via m1,b1 
        fit f2(x) "$1" using 4:5 i 2 via m2,b2 
        fit f3(x) "$1" using 4:5 i 3 via m3,b3 

        set output "together_1.pdf"

        set xtics 2
        # set xrange[-26:-2]
        # set yrange[0:8]
        set title "Volume and vertical position of cells"
        set xlabel 'z coordinate of seed (cm)'
        set ylabel 'Volume of Voronoi cell (cm^3)'

        c0 = "#0072BD"
        c1 = "#D95319"
        c2 = "#EDB120"
        c3 = "#7E2F8E"
        
        plot \
            "$OUT" i 0 using 4:5 with points lw 0.02 lt 0 lc rgb c0 notitle, \
            "$OUT" i 1 using 4:5 with points lw 0.1 lt 0 lc rgb c1 notitle, \
            "$OUT" i 2 using 4:5 with points lw 0.1 lt 0 lc rgb c2 notitle, \
            "$OUT" i 3 using 4:5 with points lw 0.1 lt 0 lc rgb c3 notitle, \
            [Z0_min:Z0_max] f0(x) with lines lw 4 lc rgb c0 title "Adult", \
            [Z0_min:Z0_max] f0(x) with lines lw 4 dt 2 lc rgb 'black' notitle, \
            [Z1_min:Z1_max] f1(x) with lines lw 4 lc rgb c1 title "13 year old", \
            [Z1_min:Z1_max] f1(x) with lines lw 4 dt 2 lc rgb 'black' notitle, \
            [Z2_min:Z2_max] f2(x) with lines lw 4 lc rgb c2 title "8 year old", \
            [Z2_min:Z2_max] f2(x) with lines lw 4 dt 2 lc rgb 'black' notitle, \
            [Z3_min:Z3_max] f3(x) with lines lw 4 lc rgb c3 title "3 year old", \
            [Z3_min:Z3_max] f3(x) with lines lw 4 dt 2 lc rgb 'black' notitle

EOF
