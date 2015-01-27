close all; clear all; clc;
load('./mat_files/hw2_kArrays.mat')

scheme = 'jet';

Ts = 1e-6;

unit10um = zeros(round(10e-6/Ts)+1,round(10e-6/Ts)+1);
unit10um(1,1) = 1;

gridsize = 400;
tmp = repmat(unit10um,[100 100]);
center_dim = round(size(tmp)/2)
dim2 = center_dim(1)+round(gridsize/2)
dim1 = center_dim(1)-round(gridsize/2)
NFFT = 2^(nextpow2(dim2-dim1))+1

im_blood = tmp(dim1:dim2, dim1:dim2);
fft_blood = fftshift(abs(fft2(im_blood,NFFT,NFFT)));

k_axes = (pi/(Ts*1000))*linspace(-1,1,NFFT);
x_axes = Ts*gridsize/2*linspace(-1,1,gridsize);

figure;
subplot(121)
imagesc(x_axes,x_axes,im_blood); axis image; colormap(scheme); 
subplot(122)
imagesc(k_axes,k_axes,fft_blood,[0 1000]); axis image; colormap(scheme);
% xlim([min(k_x{1}) max(k_x{1})])
% ylim([min(k_y{1}) max(k_y{1})])

tmpxi = find(k_x{1}>=-30 & k_x{1}<=30);
tmpyi = find(k_y{1}>=-30 & k_y{1}<=30);

[Xq,Yq] = meshgrid(k_x{1}(tmpxi),k_y{1}(tmpyi));
fft_blood_itp = interp2(k_axes,k_axes,fft_blood,Xq,Yq,'cubic');

figure
imagesc(k_x{1}(tmpxi),k_y{1}(tmpyi),fft_blood_itp,[0 1000]); axis image; colormap(scheme); 
% xlim([-30 30])
% ylim([-30 30])
% colorbar

tmp_rot = imrotate(tmp,22.5);
im_blood_rot = tmp_rot(dim1:dim2, dim1:dim2);
fft_blood_rot = fftshift(abs(fft2(im_blood_rot,NFFT,NFFT)));

figure;
subplot(121)
imagesc(x_axes,x_axes,im_blood_rot); axis image; colormap(scheme); 
subplot(122)
imagesc(k_axes,k_axes,fft_blood_rot,[0 1000]); axis image; colormap(scheme); 
% xlim([min(k_x{1}) max(k_x{1})])
% ylim([min(k_y{1}) max(k_y{1})])

fft_blood_rot_itp = interp2(k_axes,k_axes,fft_blood_rot,Xq,Yq,'cubic');

figure
imagesc(k_x{1}(tmpxi),k_y{1}(tmpyi),fft_blood_rot_itp,[0 1000]); axis image; colormap(scheme);

%%

% clear all; clc;
% load('hw2_kArrays.mat')

scheme = 'jet';

Ts = 1e-4;

unit1mm = zeros(round(1e-3/Ts)+1,round(1e-3/Ts)+1);
unit1mm(1,1) = 1;

gridsize = 300;
tmp = repmat(unit1mm,[100 100]);
center_dim = round(size(tmp)/2)
dim2 = center_dim(1)+round(gridsize/2);
dim1 = center_dim(1)-round(gridsize/2);
NFFT = 2^(nextpow2(dim2-dim1))+1

im_port = tmp(dim1:dim2, dim1:dim2);
fft_port = fftshift(abs(fft2(im_port,NFFT,NFFT)));

k_axes = (pi/(Ts*1000))*linspace(-1,1,NFFT);
x_axes = Ts*gridsize/2*linspace(-1,1,gridsize);

figure;
subplot(121)
imagesc(x_axes,x_axes,im_port); axis image; colormap(scheme);
subplot(122)
imagesc(k_axes,k_axes,fft_port,[0 900]); axis image; colormap(scheme); 
% xlim([-30 30])
% ylim([-30 30])
% xlim([min(k_x{1}) max(k_x{1})])
% ylim([min(k_y{1}) max(k_y{1})])

tmpxi = find(k_x{1}>=-30 & k_x{1}<=30);
tmpyi = find(k_y{1}>=-30 & k_y{1}<=30);

[Xq,Yq] = meshgrid(k_x{1}(tmpxi),k_y{1}(tmpyi));
fft_port_itp = interp2(k_axes,k_axes,fft_port,Xq,Yq,'cubic');
% size(fft_port)
% size(fft_port_itp)
figure
imagesc(k_x{1}(tmpxi),k_y{1}(tmpyi),fft_port_itp,[0 1000]); axis image; colormap(scheme);

% %%
% clear all; clc
% load('hw2_kArrays.mat')

scheme = 'jet';

Ts = 1e-5;
gridsize = 300;
tmp = rand(gridsize);
NFFT = 2^(nextpow2(gridsize))+1;

x_axes = Ts*gridsize/2*linspace(-1,1,gridsize);
k_axes = (pi/(Ts*1000))*linspace(-1,1,NFFT);

thresh = 0.9;
tmp(find(tmp<thresh)) = 0;
tmp(find(tmp>=thresh)) = 1;

im_rnd = tmp;
fft_rnd = fftshift(abs(fft2(im_rnd,NFFT,NFFT)));

figure
subplot(121)
imagesc(x_axes,x_axes,im_rnd); axis image; colormap(scheme);
subplot(122)
imagesc(k_axes,k_axes,fft_rnd,[0 900]); axis image; colormap(scheme); 

tmpxi = find(k_x{1}>=-30 & k_x{1}<=30);
tmpyi = find(k_y{1}>=-30 & k_y{1}<=30);

[Xq,Yq] = meshgrid(k_x{1}(tmpxi),k_y{1}(tmpyi));
fft_rnd_itp = interp2(k_axes,k_axes,fft_rnd,Xq,Yq,'cubic');

figure
imagesc(k_x{1}(tmpxi),k_y{1}(tmpyi),fft_rnd_itp,[0 1000]); axis image; colormap(scheme);

save('hw2_targets.mat','tmpxi','tmpyi','fft_rnd_itp','fft_port_itp','fft_blood_rot_itp','fft_blood_itp')