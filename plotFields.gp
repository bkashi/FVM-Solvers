set view map
set size ratio -1

set xtics axis format " "
set ytics axis format " "
unset title
splot 'Ux' matrix with image
pause -1
