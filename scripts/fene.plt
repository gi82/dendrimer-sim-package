set terminal wxt ctrl persist enhanced title "FENE"


set title   "FENE potential"
set xlabel   "r"
set ylabel   "U"
set xrange   [0:5]
set yrange   [-5:20]


FENE(x)  = -K * R **2* log ( 1 - ((x-x0)/R)**2 )

K   = 40
R   = 0.375
x0  = 1.875

plot \
  FENE(x) with lines lc rgb "#FF0000" title "FENE"    
