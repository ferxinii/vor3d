set terminal pngcairo enhanced font 'Arial,18' size 1080,1080 enhanced
set output "sphere_vol_1.png"

binwidth = 0.5;
bin(x,width)=width*floor(x/width) + width/2.0

set title "Distribution of cell volumes"
set xlabel 'Volume (??)'
set ylabel 'Count'

set xtics 2
set grid
set boxwidth binwidth
set style fill solid border -1 
set xrange[0:60]
plot 'test_volumes.txt' using (bin($4,binwidth)):(1.0) smooth freq with boxes lc rgb "#0072bd" notitle


set xrange[-26:-2]
set yrange[0:10]

set output "sphere_vol_2.png"
plot 'test_volumes.txt' using 3:4 with points notitle

