set terminal wxt
set title "Cp vs. x/c"
set xlabel "x/c"
set ylabel "Cp"
set datafile separator ','

stats 'VortexPanel.csv' using 2 nooutput
YMAX = STATS_max + 0.5
YMIN = STATS_min - 0.5

stats 'VortexPanel.csv' using 1 nooutput
XMAX = STATS_max + 0.05
XMIN = STATS_min - 0.05

set yrange [YMAX:YMIN]
set xrange [XMIN:XMAX]

plot 'VortexPanel.csv' using 1:2 with lines linewidth 3 notitle
