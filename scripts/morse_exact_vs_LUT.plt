set terminal wxt ctrl persist enhanced title "test"


set title   "Morse potential"
set xlabel   "r"
set grid
set ylabel   "U(r)"
set xrange   [-1:7]
set yrange   [-100:1000]

MORSE(x) =  eps*((exp(-a*(x-d))	-1)**2-1)

eps = 0.714
a   = 6.40
d   = 1.00

plot MORSE(x) with lines lc rgb "#00FFF0" title "Morse function"\
,'morse.dat' us 1:2 with lines title "Morse LUT"\
,'morse.dat' us 1:3 with lines title "Morse exact calc"
