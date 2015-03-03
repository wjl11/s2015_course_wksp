close all; clc; clear all;
im = double(imread('cameraman.tif'));

f = fspecial('gaussian',[10 10],2.8); % weighted average of kernel
imSmooth = imfilter(im,f);

figure
subplot(121)
imagesc(im), colormap gray; axis image
subplot(122)
imagesc(imSmooth), colormap gray; axis image


edge_f = [1,0,-1]; % lateral edge detection (will zero out when there are no edges)
imEdgeH = imfilter(im,edge_f);
figure
subplot(141)
imagesc(im), colormap gray; axis image

subplot(142)
imagesc(abs(imEdgeH)), colormap gray; axis image

imEdgeV = imfilter(im, edge_f');

subplot(143)
imagesc(abs(imEdgeV)), colormap gray; axis image

imEdgeAll = abs(imEdgeV)+abs(imEdgeH);

subplot(144)
imagesc(imEdgeAll), colormap gray; axis image

[gx,gy] = gradient(im);
figure;
subplot(122),imagesc(gx)
subplot(121),imagesc(gy)

mag = sqrt(gx.^2+gy.^2);
ang = atan2(gy,gx);

figure
subplot(121)
imagesc(mag)
subplot(122)
imagesc(ang)

% hue, value, saturation
H = (ang+pi)/(2*pi);
V = mag/max(mag(:));
S = ones(size(mag));
RGB = hsv2rgb(cat(3,H,S,V));
figure
imagesc(RGB)

imEdge = edge(im);
figure
imagesc(imEdge);

[H,T,R] = hough(imEdge);
figure
imagesc(T,R,H); colormap hot

% non-maximum supression
nmCam = nonmax(mag,ang); % NON MATLAB FUNCTION
figure
imagesc(nmCam)

% hysteresis thresholding 

