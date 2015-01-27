clear all; close all; clc;
scheme = 'gray';
load('./mat_files/hw2_psf_kspace_v2.mat')
load('./mat_files/hw2_targets_v3.mat')

x_lab = 'Lateral frequency (mm^-^1)';
y_lab = 'Axial frequency (mm^-^1)';
plot_lab = {'Focused Tx','PW Tx'};

psf_lab = {'Focused Tx PSF (k-space)','PW Tx PSF (k-space)'}

k_x{1} = k_x{1};
k_y{1} = k_y{1};

for i = 1:2
fft_psf{i} = kspace{i}(tmpyi,tmpxi);
fft_im_blood{i} = fft_psf{i}.*fft_blood_itp;
fft_im_blood_rot{i} = fft_psf{i}.*fft_blood_rot_itp;
fft_im_port{i} = fft_psf{i}.*fft_port_itp;
fft_im_rnd{i} = fft_psf{i}.*fft_rnd_itp;

figure(1)
subplot(1,2,i)
imagesc(k_x{1}(tmpxi),k_y{1}(tmpyi),fft_psf{i}); colormap(scheme); axis image
xlabel(x_lab,'Fontsize',14)
ylabel(y_lab,'Fontsize',14)
title(psf_lab{i},'Fontsize',14)
    set(gcf,'position',[100 100 1000 400])
    set(gca,'YDir','normal');


figure(2)
subplot(1,2,i)
imagesc(k_x{1}(tmpxi),k_y{1}(tmpyi),fft_im_blood{i}); colormap(scheme); axis image
xlabel(x_lab,'Fontsize',14)
ylabel(y_lab,'Fontsize',14)
title(['10 \mum grid (' plot_lab{i} ')'],'Fontsize',14)
    set(gcf,'position',[100 100 1000 400])
    set(gca,'YDir','normal');


figure(3)
subplot(1,2,i)
imagesc(k_x{1}(tmpxi),k_y{1}(tmpyi),fft_im_blood_rot{i}); colormap(scheme); axis image
xlabel(x_lab,'Fontsize',14)
ylabel(y_lab,'Fontsize',14)
title(['Rotated 10 \mum grid (' plot_lab{i} ')'],'Fontsize',14)
    set(gcf,'position',[100 100 1000 400])
    set(gca,'YDir','normal');


figure(4)
subplot(1,2,i)
imagesc(k_x{1}(tmpxi),k_y{1}(tmpyi),fft_im_port{i}); colormap(scheme); axis image
xlabel(x_lab,'Fontsize',14)
ylabel(y_lab,'Fontsize',14)
title(['1 mm grid (' plot_lab{i} ')'],'Fontsize',14)
    set(gcf,'position',[100 100 1000 400])
    set(gca,'YDir','normal');


figure(5)
subplot(1,2,i)
imagesc(k_x{1}(tmpxi),k_y{1}(tmpyi),fft_im_rnd{i}); colormap(scheme); axis image;
xlabel(x_lab,'Fontsize',14)
ylabel(y_lab,'Fontsize',14)
title(['Random grid (' plot_lab{i} ')'],'Fontsize',14)
    set(gcf,'position',[100 100 1000 400])
    set(gca,'YDir','normal');

end

figure(6)
imagesc(ifftshift(abs(ifft(fft_im_port{1})),1))

% figure(1)
%     print -dpng -r300 ./imkspace_psf.png
% figure(2)
%     print -dpng -r300 ./imkspace_10um.png
% figure(3)
%     print -dpng -r300 ./imkspace_10um_rot.png
% figure(4)
%     print -dpng -r300 ./imkspace_1mm.png
% figure(5)
%     print -dpng -r300 ./imkspace_rnd.png