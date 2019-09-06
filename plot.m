M=dlmread("output.txt"," ");
h=imagesc(flipud(M'));
axis equal
axis on

colorbar("South")
colormap(jet)
%set(h,'Position',[60 300 1200 500])
pause

