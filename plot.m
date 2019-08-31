M=dlmread("output.txt"," ");
surf(M','edgecolor','none')
axis equal
axis off
colorbar
view(0,90)
colormap(jet)
pause
