set terminal wxt ctrl persist enhanced title "Morse Potentials" 0

set key
set title   "Morse potential"
set xlabel   "r"
set ylabel   "U(r)"
#set xtics    0.2	
set xrange   [0:4]
set yrange   [-5:5]
set grid

MORSE1(x) =  eps*((exp(-a*(x-d))	-1)**2-1)
MORSE2(x) =  eps*((exp(-a*(x-d))	-1)**2)

set multiplot
eps = 0.714
a   = 6.4
d   = 1.00
plot MORSE1(x) with lines lc rgb "#00FFF0" title "CC"

eps = 0.714
a   = 10.4
d   = 1.00
plot MORSE1(x) with lines lc rgb "#00GGG0" title"CS"

unset multiplot


