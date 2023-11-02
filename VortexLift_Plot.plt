set terminal wxt
set title "Cl vs. Theta ({/Symbol q})"
set xlabel "Theta ({/Symbol q})"
set ylabel "Cl"
set datafile separator ','

stats 'Lift.csv' using 2 nooutput
YMAX = STATS_max + 2
YMIN = STATS_min - 2

stats 'Lift.csv' using 1 nooutput
XMAX = STATS_max + 5
XMIN = STATS_min - 5

set yrange [YMIN:YMAX]
set xrange [XMIN:XMAX]

plot 'Lift.csv' using 1:2 with linespoints linewidth 3 notitle, '' using 1:2:(sprintf("(%g,%g)", $1, $2)) with labels offset 0,1 notitle
