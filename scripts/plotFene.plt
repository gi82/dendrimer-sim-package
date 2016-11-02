set terminal wxt ctrl persist enhanced title "FENE potentials"

set title   "FENE potential"
set xlabel   "r"
#set xtics 0.5
set ylabel   "U"
#set xrange   [-10:10]
#set yrange   [0:20]


FENE1(x)  = -K1 * R1 **2* log ( 1 - ((x-x01)/R1)**2 )
#FENE2(x)  = -K2 * R2 **2* log ( 1 - ((x-x02)/R2)**2 )
#FENE3(x)  = -K3 * R3 **2* log ( 1 - ((x-x03)/R3)**2 )
#CC params
K1   = 40
x01  = 2
R1   = 1

#CS params
K2   = 30
x02  = 3.750
R2   = 0.750

#BB params
K3   = 40
R3   = 2.8125
x03  = 0.5625

plot\
      FENE1(x) with lines lc rgb "#00AAA0" title"CC"
#     ,FENE2(x) with lines lc rgb "#00AGG0" title"CS"\
#     ,FENE3(x) with lines lc rgb "#00AGA0" title"BB"
