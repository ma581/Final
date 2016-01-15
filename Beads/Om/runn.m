% Used to show the heat map generated with Om's code

close all;
clear;
filename = 'C:\Users\Manoj\Dropbox\SainsburyLab\Beads\Om\2015_05_03_beads_3';
rad = 100;
IlluminationFlatness(filename,rad);

close all
load('beadsCalib100.mat')
h = surf(beadsCalib);
view(2);
colorbar;
title('Mean Intensity')           
set(h,'LineStyle','none')


%  set(gca,'clim',[minval maxval])