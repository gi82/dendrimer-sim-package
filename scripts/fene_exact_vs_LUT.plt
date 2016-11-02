set terminal wxt ctrl persist enhanced title "test"


set title   "Fene potential"
set xlabel   "r"
set grid
set ylabel   "Fene(r)"
set xrange   [-1:7]
set yrange   [-1:50]

FENE1(x)  = -K1 * R1 **2* log ( 1 - ((x-x01)/R1)**2 )
#CC params
K1   = 40
x01  = 1.875
R1   = 0.375

#CS params
K2   = 30
x02  = 3.750
R2   = 0.750

#BB params
K3   = 40
R3   = 2.8125
x03  = 0.5625


plot FENE1(x) with lines lc rgb "#00FFF0" title "Morse function"\
,'fene.dat' us 1:2 with lines title "Morse LUT"\
#,'fene.dat' us 1:3 with lines title "Morse exact calc"
