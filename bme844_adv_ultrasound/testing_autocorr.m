clear all; clc; 
close all;

load('./mat_files/hw5_data_full.mat');

dz = (1/fs)*(1540/2);
depth = foc_z;

rf = rf{10,10};
env = env{10,10};

figure; imagesc(rf);
figure; imagesc(env);

rfk = rf(400:600,:);
envk = env(400:600,:);

figure
imagesc(angle(envk)); colorbar

rfcorr = normxcorr2(rfk./max(rfk(:)),rfk./max(rfk(:)));
envcorr = normxcorr2(envk./max(envk(:)),envk./max(envk(:)));

figure; imagesc(rfcorr);
figure; imagesc(envcorr);

figure; hold on
plot(rfcorr(ceil(size(rfcorr,1)/2),:))
plot(envcorr(ceil(size(envcorr,1)/2),:),'r')
hold off
legend('rf','detected')