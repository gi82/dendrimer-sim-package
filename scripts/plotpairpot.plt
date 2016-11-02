set encoding iso_8859_1
set terminal postscript fontfile "/usr/share/texmf/fonts/type1/public/lm/lmmi10.pfb"
set terminal postscript fontfile "/usr/share/texmf/fonts/type1/public/lm/lmr10.pfb"
set terminal postscript fontfile "/usr/share/texmf/fonts/type1/public/lm/lmsy10.pfb"

rm = "LMRoman10-Regular"       ; rm ( text ) = "{/".rm." ".text."}"
mi = "LMMathItalic10-Regular"  ; mi ( text ) = "{/".mi." ".text."}"
sy = "LMMathSymbols10-Regular" ; sy ( text ) = "{/".sy." ".text."}"

unset key

set terminal postscript color enhanced font rm.",20" eps
set output "pairpotD5D7.eps"

set border lw 2

set xtics (\
   "0"  0, "" 0.5 1, \
   "1"  1, "" 1.5 1, \
   "2"  2, "" 2.5 1, \
   "3"  3, "" 3.5 1, \
   "4"  4, "" 4.5 1 \
)

set ytics (\
   "-0.5"  -0.5,"" 0.25 1,\
   "0"  0,"" 0.25 1,\
   "0.5"  0.5,"" 0.75 1,\
   "1"  1, "" 1.25 1,\
   "1.5"  1.5,"" 1.75 1,\
   "2"  2, "" 2.25 1,\
   "2.5"  2.5, "" 2.75 1,\
   "3"  3\
)


set xlabel   " (r/R_{g})" #offset 1.9, 0.0	
set ylabel   mi("\321\315")." (r)"  #offset 1.9, 0.0	

set xrange   [0:4.5]
set yrange   [-0.5:3]

set format x ""
set format y ""
f(x)=a
f1(x)=b
fc(x)=0
fit f(x) "<awk '{if (NR>85){print $4,$7}}' potD5.hs.100.dat" via a
fit f1(x) "<awk '{if (NR>85){print $4,$7}}' potD7.hs.100.dat" via b


plot \
  "potD5.hs.100.dat" us 4:($7-a) with points lw 4 \
  ,"potD7.hs.100.dat" us 4:($7-b) with points lw 4 \
   ,fc(x)
#   ,f1(x)
