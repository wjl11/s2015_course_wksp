clear all; clc;

load('hw2_foc_psf_v2.mat');

figure
for i = 1:2
    subplot(1,2,i)
    imagesc(x*100,z{i}*100,rf{i}); axis image; colormap jet
end

for i = 1:2
   
   N_y = size(rf{i},1) % number of image rows (y aspect)
   N_x = size(rf{i},2) % number of image columns
   NFFT_x = 2^nextpow2(N_x)+1
   NFFT_y = 2^nextpow2(N_y)+1
   
   kspace{i} = fftshift(abs(fft2(rf{i},NFFT_y,NFFT_x)));
   kspace_nopad{i} = fftshift(abs(fft2(rf{i})));

   tmp = diff(z{i});
   T_y(i) = 1000*tmp(1)
   tmp = diff(x);
   T_x(i) = 1000*tmp(1)
%    k_tmp_y = 2*pi*(0:(N_y-1))/N_y;
%    k_tmp_x = 2*pi*(0:(N_x-1))/N_x;
%    k_y2{i} = unwrap(fftshift(k_tmp_y)-2*pi)./(T_y(i));
%    k_x2{i} = unwrap(fftshift(k_tmp_x)-2*pi)./(T_x(i));
   
   k_y{i} = (1/(2*T_y(i)))*linspace(-1,1,NFFT_y);
   k_x{i} = (1/(2*T_x(i)))*linspace(-1,1,NFFT_x);
   
%    k_y_nopad{i} = (pi/T_y(i))*linspace(-1,1,N_y);
%    k_x_nopad{i} = (pi/T_x(i))*linspace(-1,1,N_x);
   
   clear N_y N_x k_tmp_y k_tmp_x tmp
end

figure
for i = 1:2
    subplot(1,2,i)
    imagesc(k_x{i},k_y{i},kspace{i}); colormap jet; axis image; 
%     colorbar
    xlim([-5 5])
    ylim([-10 10])
    xlabel('Lateral frequency (mm^-^1)')
    ylabel('Axial frequency (mm^-^1)')
    max(k_y{i})
    max(k_x{i})
end

% figure
% for i = 1:2
%     subplot(1,2,i)
%     imagesc(k_x_nopad{i},k_y_nopad{i},kspace_nopad{i}); colormap jet; axis image; colorbar
% %     xlim([-30 30])
% %     ylim([-30 30])
%     xlabel('k_l_a_t (rads/mm)')
%     ylabel('k_a_x (rads/mm)')
% end

% save('hw2_psf_kspace_v2.mat','k_y','k_x','kspace')





