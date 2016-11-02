set terminal wxt ctrl persist enhanced title "test"


set title   "Morse potential"
set xlabel   "r"
set grid
set ylabel   "U(r)"
#set xrange   [-1:7]
#set yrange   [-100:1000]
plot 'feneCS.dat' us 1:2 with lines title "Morse LUT"
