clear all; close all; clc

figure
subplot(121)
load bmode_1.0.mat
imagesc(bws*1e3,deps*1e3,env,[-40 0]), colormap gray
xlabel('y (mm)'), ylabel('z (mm)')
axis equal, axis tight
title('1 MHz')
% cbar=colorbar; title(cbar,'dB')

subplot(122)
load bmode_1.5.mat
imagesc(bws*1e3,deps*1e3,env,[-40 0]), colormap gray
xlabel('y (mm)'), ylabel('z (mm)')
axis equal, axis tight
title('1.5 MHz')
% cbar=colorbar; title(cbar,'dB')

print -djpeg fig2b_WillieLong.jpg