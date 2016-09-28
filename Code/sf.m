function filename = sf(str)
%save figure to Img folder

filename = sprintf(['./images/',str,'.png']);
saveas(gcf, filename,'png');
fprintf(['Saved Results to ' filename '\n']);