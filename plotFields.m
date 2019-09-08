Ux=dlmread("Ux"," ");
Uy=dlmread("Uy"," ");
M=dlmread("T"," ");

h=figure;
axis equal

imagesc(flipud(M'));
axis equal

hold on
quiver( Ux, -Uy, 'k' );
axis equal

set(h,'Position',[60 200 1800 800])

axis off
%colorbar("South")
colormap(jet)

pause

