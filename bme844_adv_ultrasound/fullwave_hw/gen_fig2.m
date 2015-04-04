clear all; close all; clc

figure
subplot(121)
load PSF_homogenous_-13.mat
imagesc(lats*1e3,deps*1e3,env,[-45 0]), colormap gray
xlabel('mm'),ylabel('mm')
axis equal, axis tight
title(sprintf('R = %.3f', R))
cbar=colorbar; title(cbar,'dB')

subplot(122)
load PSF_homogenous_500.mat
imagesc(lats*1e3,deps*1e3,env,[-45 0]), colormap gray
xlabel('mm'),ylabel('mm')
axis equal, axis tight
title(sprintf('R = %.1f', R))
cbar=colorbar; title(cbar,'dB')
print -djpeg fig2a_1_WillieLong.jpg

figure
subplot(121)
load PSF_tissue_-13.mat
imagesc(lats*1e3,deps*1e3,env,[-45 0]), colormap gray
xlabel('mm'),ylabel('mm')
axis equal, axis tight
title(sprintf('R = %.3f', R))
cbar=colorbar; title(cbar,'dB')

subplot(122)
load PSF_tissue_500.mat
imagesc(lats*1e3,deps*1e3,env,[-45 0]), colormap gray
xlabel('mm'),ylabel('mm')
axis equal, axis tight
title(sprintf('R = %.1f', R))
cbar=colorbar; title(cbar,'dB')
print -djpeg fig2a_2_WillieLong.jpg