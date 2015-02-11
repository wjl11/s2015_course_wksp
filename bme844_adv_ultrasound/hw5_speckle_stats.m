clear all; clc; 
close all;

load('./mat_files/hw5_data_full.mat');

dz = (1/fs)*(1540/2);
depth = foc_z;

% ax_res = 0.75*lambda;

for n = 1:length(scatPerResCell)
    for i = 1:length(rndSeed)
        n_r{n,i} = find(r{n,i}>zposlim(1)+.8*lambda & r{n,i}<zposlim(2)-.8*lambda);
        
        n_snr{n,i} = find(r{n,i}>zposlim(1)+2.5*lambda & r{n,i}<zposlim(2)-2.5*lambda);
        
        trim_env = env{n,i}(n_r{n,i},:);
        trim_rf = rf{n,i}(n_r{n,i},:);
        env_snr = env{n,i}(n_snr{n,i},:);
        SNR(n,i) = mean(env_snr(:))/std(env_snr(:));
           
        center_i = size(trim_env,2);
        center_j = size(trim_env,1);
        
        n_ave = 3;
        corr_i = (center_i-floor(n_ave/2)):(center_i+floor(n_ave/2));
        corr_j = (center_j-floor(n_ave/2)):(center_j+floor(n_ave/2));
        
        corr_env{n,i} = normxcorr2(trim_env/max(trim_env(:)), trim_env/max(trim_env(:)));
        corr_rf{n,i} = normxcorr2(trim_rf/max(trim_rf(:)), trim_rf/max(trim_rf(:)));
        
        dz = (1/fs)*(1540/2);
        z = (0:center_j-1).*dz;
        x = linspace(0,max(xposlim),center_i);
        
        thresh = 0.45;
        for nn = 1:n_ave
            tmp_lat_rf = abs(hilbert(corr_rf{n,i}(corr_j(nn),:)));
            tmp_lat_rf = tmp_lat_rf/max(tmp_lat_rf);
            tmp_lat_env = corr_env{n,i}(corr_j(nn),:);
            tmp_lat_env = tmp_lat_env/max(tmp_lat_env);

            tmp_ax_rf = abs(hilbert(corr_rf{n,i}(:,corr_i(nn))));
            tmp_ax_rf = tmp_ax_rf/max(tmp_ax_rf);
            tmp_ax_env = corr_env{n,i}(:,corr_i(nn));
            tmp_ax_env = tmp_ax_env/max(tmp_ax_env);

            raw_lat_rf(nn) = 2*x(find(tmp_lat_rf(center_i:end) <= thresh,1));
            raw_lat_env(nn) = 2*x(find(tmp_lat_env(center_i:end) <= thresh,1));
            raw_ax_rf(nn) = 2*z(find(tmp_ax_rf(center_j:end) <= thresh,1));
            raw_ax_env(nn) = 2*z(find(tmp_ax_env(center_j:end) <= thresh,1));
            
%             keyboard

            clear tmp_lat_rf tmp_lat_env tmp_ax_rf tmp_ax_env
            
            
        end
        lat_res_rf(n,i) = mean(raw_lat_rf);
        lat_res_env(n,i) = mean(raw_lat_env);
        ax_res_rf(n,i) = mean(raw_ax_rf);
        ax_res_env(n,i) = mean(raw_ax_env);
        clear raw_lat_rf raw_lat_env raw_ax_rf raw_ax_env

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

