set terminal png enh size 1200,600 enhanced font "Times New Roman,12"
set output 'solution.png'
set nokey
set pointsize 0.25
set xlabel "x"
set ylabel "y"
set view 60,30,1
splot \
 'solution.dat' using 5:6:7 title "" with points pt 7
