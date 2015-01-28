close all; clear all; clc;
load('./mat_files/hw2_psf_kspace_v2.mat')

scheme = 'gray';
x_labx = 'Lateral Position (mm)';
y_labx = 'Axial Position (mm)';
x_labk = 'Lateral frequency (mm^-^1)';
y_labk = 'Axial frequency (mm^-^1)';

% Ts = 1e-5; % use for interp data
Ts = 5e-7; % use for visualization data`

unit10um = zeros(round(10e-6/Ts)+1,round(10e-6/Ts)+1);
unit10um(1,1) = 1;

gridsize = 300;
tmp = repmat(unit10um,[400 400]);
center_dim = round(size(tmp)/2)
dim2 = center_dim(1)+round(gridsize/2)
dim1 = center_dim(1)-round(gridsize/2)
NFFT = dim2-dim1+1
% 2^(nextpow2(dim2-dim1)+1)+1

im_blood = tmp(dim1:dim2, dim1:dim2);
fft_blood = fftshift(abs(fft2(im_blood,NFFT,NFFT)));

k_axes = (1/(2*Ts*1000))*linspace(-1,1,NFFT);
x_axes = Ts*gridsize/2*linspace(-1,1,gridsize);

figure;
subplot(121)
imagesc(1000*x_axes,1000*x_axes,im_blood); axis image; colormap(scheme); 
xlabel(x_labx,'Fontsize',14)
ylabel(y_labx,'Fontsize',14)
title('10 \mum point grid','Fontsize',14)
subplot(122)
imagesc(k_axes,k_axes,fft_blood); axis image; colormap(scheme);
xlabel(x_labk,'Fontsize',14)
ylabel(y_labk,'Fontsize',14)
title('10 \mum point grid (k-space)','Fontsize',14)
set(gcf,'position',[100 100 1000 400])

% print -dpng -r300 ./kspace_10um.png
% xlim([-30 30])
% ylim([-30 30])

tmpxi = find(k_x{1}>=-5 & k_x{1}<=5);
tmpyi = find(k_y{1}>=-10 & k_y{1}<=10);

[Xq,Yq] = meshgrid(k_x{1}(tmpxi),k_y{1}(tmpyi));
fft_blood_itp = interp2(k_axes,k_axes,fft_blood,Xq,Yq,'cubic');

figure
imagesc(k_x{1}(tmpxi),k_y{1}(tmpyi),fft_blood_itp); axis image; colormap(scheme); 

tmp_rot = imrotate(tmp,22.5);
im_blood_rot = tmp_rot(dim1:dim2, dim1:dim2);
fft_blood_rot = fftshift(abs(fft2(im_blood_rot,NFFT,NFFT)));

figure;
subplot(121)
imagesc(1000*x_axes,1000*x_axes,im_blood_rot); axis image; colormap(scheme); 
xlabel(x_labx,'Fontsize',14)
ylabel(y_labx,'Fontsize',14)
title('Rotated 10 \mum point grid','Fontsize',14)
subplot(122)
imagesc(k_axes,k_axes,fft_blood_rot); axis image; colormap(scheme); 
xlabel(x_labk,'Fontsize',14)
ylabel(y_labk,'Fontsize',14)
title('Rotated 10 \mum point grid (k-space)','Fontsize',14)
set(gcf,'position',[100 100 1000 400])

% print -dpng -r300 ./kspace_10um_rot.png
% xlim([min(k_x{1}) max(k_x{1})])
% ylim([min(k_y{1}) max(k_y{1})])
% xlim([-30 30])
% ylim([-30 30])

fft_blood_rot_itp = interp2(k_axes,k_axes,fft_blood_rot,Xq,Yq,'cubic');

figure
imagesc(k_x{1}(tmpxi),k_y{1}(tmpyi),fft_blood_rot_itp); axis image; colormap(scheme);

%%

% clear all; clc;
% load('hw2_kArrays.mat')

scheme = 'gray';

Ts = 2e-5;

unit1mm = zeros(round(1e-3/Ts)+1,round(1e-3/Ts)+1);
unit1mm(1,1) = 1;

gridsize = 300;
tmp = repmat(unit1mm,[300 300]);
center_dim = round(size(tmp)/2)
dim2 = center_dim(1)+round(gridsize/2);
dim1 = center_dim(1)-round(gridsize/2);
NFFT = dim2-dim1+1
% 2^(nextpow2(dim2-dim1))+1

im_port = tmp(dim1:dim2, dim1:dim2);
fft_port = fftshift(abs(fft2(im_port,NFFT,NFFT)));

k_axes = (1/(2*Ts*1000))*linspace(-1,1,NFFT);
x_axes = Ts*gridsize/2*linspace(-1,1,gridsize);

figure;
subplot(121)
imagesc(1000*x_axes,1000*x_axes,im_port); axis image; colormap(scheme);
xlabel(x_labx,'Fontsize',14)
ylabel(y_labx,'Fontsize',14)
title('1 mm point grid','Fontsize',14)
subplot(122) 
imagesc(k_axes,k_axes,fft_port); axis image; colormap(scheme); 
xlabel(x_labk,'Fontsize',14)
ylabel(y_labk,'Fontsize',14)
title('1 mm point grid (k-space)','Fontsize',14)
set(gcf,'position',[100 100 1000 400])

% print -dpng -r300 ./kspace_1mm.png

% xlim([-30 30])
% ylim([-30 30])
% xlim([min(k_x{1}) max(k_x{1})])
% ylim([min(k_y{1}) max(k_y{1})])

tmpxi = find(k_x{1}>=-5 & k_x{1}<=5);
tmpyi = find(k_y{1}>=-10 & k_y{1}<=10);

[Xq,Yq] = meshgrid(k_x{1}(tmpxi),k_y{1}(tmpyi));
fft_port_itp = interp2(k_axes,k_axes,fft_port,Xq,Yq,'cubic');
% size(fft_port)
% size(fft_port_itp)

figure
imagesc(k_x{1}(tmpxi),k_y{1}(tmpyi),fft_port_itp); axis image; colormap(scheme);

%%
% clear all; clc
% load('hw2_kArrays.mat')

scheme = 'gray';

Ts = 5e-5;
gridsize = 300;
tmp = rand(gridsize);
NFFT = gridsize;
% 2^(nextpow2(gridsize))+1;

x_axes = Ts*gridsize/2*linspace(-1,1,gridsize);
k_axes = (1/(2*Ts*1000))*linspace(-1,1,NFFT);

thresh = 0.6;
tmp(find(tmp<thresh)) = 0;
tmp(find(tmp>=thresh)) = 1;

im_rnd = tmp;
fft_rnd = fftshift(abs(fft2(im_rnd,NFFT,NFFT)));

figure
subplot(121)
imagesc(1000*x_axes,1000*x_axes,im_rnd); axis image; colormap(scheme);
xlabel(x_labx,'Fontsize',14)
ylabel(y_labx,'Fontsize',14)
title('Randomly spaced point grid','Fontsize',14)

subplot(122)
imagesc(k_axes,k_axes,fft_rnd,[0 500]); axis image; colormap(scheme);
xlabel(x_labk,'Fontsize',14)
ylabel(y_labk,'Fontsize',14)
title('Randomly spaced point grid (k-space)','Fontsize',14)

set(gcf,'position',[100 100 1000 400])

% print -dpng -r300 ./kspace_rnd.png

tmpxi = find(k_x{1}>=-5 & k_x{1}<=5);
tmpyi = find(k_y{1}>=-10 & k_y{1}<=10);

[Xq,Yq] = meshgrid(k_x{1}(tmpxi),k_y{1}(tmpyi));
fft_rnd_itp = interp2(k_axes,k_axes,fft_rnd,Xq,Yq,'cubic');

figure
imagesc(k_x{1}(tmpxi),k_y{1}(tmpyi),fft_rnd_itp,[0 500]); axis image; colormap(scheme);
%%
% save('hw2_targets_v3.mat','tmpxi','tmpyi','fft_rnd_itp','fft_port_itp','fft_blood_rot_itp','fft_blood_itp')