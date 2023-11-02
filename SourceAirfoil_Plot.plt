set terminal wxt
set title "Cp vs x" 
set xlabel "x"
set ylabel "Cp"
set datafile separator ',' 

stats 'lower.csv' using 2 nooutput
YMAX = STATS_max + 0.1
YMIN = STATS_min - 0.1

stats 'lower.csv' using 1 nooutput
XMAX = STATS_max + 0.1
XMIN = STATS_min - 0.1

set yrange [YMAX:YMIN]
set xrange [XMIN:XMAX]


plot 'lower.csv' using 1:2 with lines pointsize 2 linewidth 3 linecolor rgb 'red' title "lower", \
     'upper.csv' using 1:2 with lines pointsize 2 linewidth 3 linecolor rgb 'blue' title "upper"

