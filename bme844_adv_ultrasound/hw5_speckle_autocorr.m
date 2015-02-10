clear all; clc; 
close all;

load('./mat_files/hw5_data_full.mat');

dz = (1/fs)*(1540/2);
depth = foc_z+0.0005;

for n = 1:length(scatPerResCell)
    for i = 1:length(rndSeed)
        n_r{n,i} = find(r{n,i}>zposlim(1)+2.5*ax_res & r{n,i}<zposlim(2)-2.5*ax_res);
        n_th{n,i} = find(rad2deg(th_scan)>-10 & rad2deg(th_scan)<10);
        
        if i == 1
        figure; 
        imagesc(th_scan(n_th{n,i}),r{n,i}(n_r{n,i}),env_db{n,i}(n_r{n,i},n_th{n,i}),[-40 0]); colormap gray
        end
        trim_env = env{n,i}(n_r{n,i},n_th{n,i});
        SNR(n,i) = mean(trim_env(:))/std(trim_env(:));
        
%         trim_rf = rf{n,i}(n_r{n,i},n_th{n,i});
        
%         corr_env{n,i} = xcorr2(trim_env./max(trim_env(:)));
%         corr_rf{n,i} = xcorr2(trim_rf./max(trim_rf(:)));
        
        corr_i = find(th_scan == 0);
        corr_j = round(length(r{n,i})/2);
%         find(r{n,i} >= depth,1);

        [ax_corr_env,ax_corr_lag] = xcorr(env{n,i}(n_r{n,i},corr_i),'coeff');
        ax_corr_rf = abs(hilbert(xcorr(rf{n,i}(n_r{n,i},corr_i),'coeff')));
        
        ax_idx = find(ax_corr_lag >= 1);
        
        ax_res_rf(n,i) = 2*dz*find(ax_corr_rf(ax_idx)<=0.48,1);
        ax_res_env(n,i) = 2*dz*find(ax_corr_env(ax_idx)<=0.48,1);
        
        [lat_corr_env,lat_corr_lag] = xcorr(env{n,i}(corr_j,n_th{n,i}),'coeff');
        lat_corr_rf = abs(hilbert(xcorr(rf{n,i}(corr_j,n_th{n,i}),'coeff')));
        
        lat_idx = find(lat_corr_lag >= 1);
        th_lat_rf = th_int*find(lat_corr_rf(lat_idx) <=0.48,1);
        th_lat_env = th_int*find(lat_corr_env(lat_idx) <=0.48,1);
        
        lat_res_rf(n,i) = depth*sind(th_lat_rf)/cosd(th_lat_rf/2);
        lat_res_env(n,i) = depth*sind(th_lat_env)/cosd(th_lat_env/2);
        
%         th = th_int.*lat_corr_lag(lat_idx);
%         x = depth.*sind(th)./cosd(th./2);
%         figure
%         plot(x,lat_corr_rf(lat_idx))
%         figure;
%         plot(x,lat_corr_env(lat_idx))
        
% 
%         tmp_rf = rf{n,i};
%         corr_rf = xcorr2(tmp_rf./max(tmp_rf(:)));
%         imagesc(1:size(corr_rf,2),,corr_rf)
%         figure;
%         hold on
%         plot(ax_corr_lag(ax_idx).*dz,ax_corr_rf(ax_idx))
%         plot(ax_corr_lag(ax_idx).*dz,ax_corr_rf(ax_idx))
%         hold off
%         figure;
%         plot(ax_corr_lag(ax_idx).*dz,ax_corr_env(ax_idx))
    end
end

SNR_m = mean(SNR,2);
SNR_s = std(SNR,[],2);

lat_rf_m = mean(lat_res_rf,2);
lat_rf_s = std(lat_res_rf,[],2);

lat_env_m = mean(lat_res_env,2);
lat_env_s = std(lat_res_env,[],2);

ax_rf_m = mean(ax_res_rf,2);
ax_rf_s = std(ax_res_rf,[],2);

ax_env_m = mean(ax_res_env,2);
ax_env_s = std(ax_res_env,[],2);

figure, hold on;
errorbar(scatPerResCell,SNR_m,SNR_s)
grid on, xlabel('Average Scatterers Per Resolution Cell')
ylabel('Speckle SNR (\mu/\sigma)')
hold off

figure, hold on;
errorbar(scatPerResCell,1000*lat_rf_m,1000*lat_rf_s)
errorbar(scatPerResCell,1000*lat_env_m,1000*lat_env_s,'r')
plot(scatPerResCell,1000*lat_res.*ones(length(scatPerResCell),1),'k--')
grid on, xlabel('Average Scatterers Per Resolution Cell')
ylabel('Lateral Speckle Size (mm)')
hold off

figure, hold on;
errorbar(scatPerResCell,1000*ax_rf_m,1000*ax_rf_s)
errorbar(scatPerResCell,1000*ax_env_m,1000*ax_env_s,'r')
plot(scatPerResCell,1000*ax_res.*ones(length(scatPerResCell),1),'k--')
grid on, xlabel('Average Scatterers Per Resolution Cell')
ylabel('Axial Speckle Size (mm)')
hold off