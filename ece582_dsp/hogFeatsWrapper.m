im = double(imread('cameraman.tif'));
[gx,gy] = gradient(im);
mag = sqrt(gx.^2+gy.^2);
ang = atan2(gy,gx);
hogFeats = computeHog(mag,ang);