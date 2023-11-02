set terminal wxt
set title "Cp vs Theta ({/Symbol q})" 
set xlabel "{/Symbol q}"
set ylabel "Cp"
set datafile separator ',' 

stats 'SourcePanel.csv' using 2 nooutput
YMAX = STATS_max + 0.5
YMIN = STATS_min - 0.5

stats 'SourcePanel.csv' using 1 nooutput
XMAX = STATS_max + 0.5
XMIN = STATS_min - 0.5

set yrange [-4:3]
set xrange [-0.5:7]

set xtics ("0" 0, "{/Symbol p}/4" pi/4, "{/Symbol p}/2" pi/2, "3{/Symbol p}/4" 3*pi/4, "{/Symbol p}" pi, "5{/Symbol p}/4" 5*pi/4, "3{/Symbol p}/2" 3*pi/2, "7{/Symbol p}/4" 7*pi/4, "2{/Symbol p}" 2*pi)

plot 'SourcePanel.csv' using 1:2 with points pointsize 2 linewidth 3 linecolor rgb 'red' title "Numerical Result", \
     'Cp_ideal_values.csv' using 1:2 with lines linecolor rgb 'black' title "Analytical Result"

