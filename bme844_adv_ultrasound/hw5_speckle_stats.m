clear all; clc; 
close all;

load('./mat_files/hw5_data_full.mat');

dz = (1/fs)*(1540/2);
depth = foc_z;

thresh = 0.49;
n_ave_i = 11;
n_ave_j = 47;

% ax_res = 0.75*lambda;
clear n_r
for n = 1:length(scatPerResCell)
    for i = 1:length(rndSeed)
        n_r{n,i} = find(r{n,i}>zposlim(1)+.75*lambda & r{n,i}<zposlim(2)-.75*lambda);
        
        n_snr{n,i} = find(r{n,i}>zposlim(1)+2.5*lambda & r{n,i}<zposlim(2)-2.5*lambda);
        
        trim_env = env{n,i}(n_r{n,i},:);
        trim_rf = rf{n,i}(n_r{n,i},:);
        env_snr = env{n,i}(n_snr{n,i},:);
        SNR(n,i) = mean(env_snr(:))/std(env_snr(:));
           
        center_i = size(trim_env,2);
        center_j = size(trim_env,1);
        
%         n_ave_i = 11;
%         n_ave_j = 47;

        corr_i = (center_i-floor(n_ave_i/2)):(center_i+floor(n_ave_i/2));
        corr_j = (center_j-floor(n_ave_j/2)):(center_j+floor(n_ave_j/2));
        
        corr_env{n,i} = normxcorr2(trim_env/max(trim_env(:)), trim_env/max(trim_env(:)));
        corr_rf{n,i} = normxcorr2(trim_rf/max(trim_rf(:)), trim_rf/max(trim_rf(:)));
        
        dz = (1/fs)*(1540/2);
        z = (0:center_j-1).*dz;
        x = linspace(0,max(xposlim),center_i);
        

%         for nn = 1:n_ave_j
%             tmp_lat_rf = abs(hilbert(corr_rf{n,i}(corr_j(nn),:)));
%             tmp_lat_rf = tmp_lat_rf/max(tmp_lat_rf);
%             tmp_lat_env = corr_env{n,i}(corr_j(nn),:);
%             tmp_lat_env = tmp_lat_env/max(tmp_lat_env);
% 
%             raw_lat_rf(nn) = 2*x(find(tmp_lat_rf(center_i:end) <= thresh,1));
%             raw_lat_env(nn) = 2*x(find(tmp_lat_env(center_i:end) <= thresh,1));
%             
% %             keyboard
% 
%             clear tmp_lat_rf tmp_lat_env  
%         end
%         for nn = 1:n_ave_i
%             tmp_ax_rf = abs(hilbert(corr_rf{n,i}(:,corr_i(nn))));
%             tmp_ax_rf = tmp_ax_rf/max(tmp_ax_rf);
%             tmp_ax_env = corr_env{n,i}(:,corr_i(nn));
%             tmp_ax_env = tmp_ax_env/max(tmp_ax_env);
%             
%             raw_ax_rf(nn) = 2*z(find(tmp_ax_rf(center_j:end) <= thresh,1));
%             raw_ax_env(nn) = 2*z(find(tmp_ax_env(center_j:end) <= thresh,1));
%             
%             clear tmp_ax_rf tmp_ax_env 
%         end
%         
%         lat_res_rf(n,i) = mean(raw_lat_rf);
%         lat_res_env(n,i) = mean(raw_lat_env);
%         ax_res_rf(n,i) = mean(raw_ax_rf);
%         ax_res_env(n,i) = mean(raw_ax_env);
%         clear raw_lat_rf raw_lat_env raw_ax_rf raw_ax_env
        
        for nn = 1:n_ave_j
            tmp_lat_rf(nn,:) = abs(hilbert(corr_rf{n,i}(corr_j(nn),:)));
            tmp_lat_rf(nn,:) = tmp_lat_rf(nn,:)/max(tmp_lat_rf(nn,:));
            tmp_lat_env(nn,:) = corr_env{n,i}(corr_j(nn),:);
            tmp_lat_env(nn,:) = tmp_lat_env(nn,:)/max(tmp_lat_env(nn,:));
        end
        
        raw_lat_rf = mean(tmp_lat_rf,1);
        raw_lat_env = mean(tmp_lat_env,1);
        lat_res_rf(n,i) = 2*x(find(raw_lat_rf(center_i:end) <= thresh,1));
        lat_res_env(n,i) = 2*x(find(raw_lat_env(center_i:end) <= thresh,1));
        clear tmp_lat_rf tmp_lat_env raw_lat_rf raw_lat_env
        
        for nn = 1:n_ave_i
            tmp_ax_rf(nn,:) = abs(hilbert(corr_rf{n,i}(:,corr_i(nn))));
            tmp_ax_rf(nn,:) = tmp_ax_rf(nn,:)/max(tmp_ax_rf(nn,:));
            tmp_ax_env(nn,:) = corr_env{n,i}(:,corr_i(nn));
            tmp_ax_env(nn,:) = tmp_ax_env(nn,:)/max(tmp_ax_env(nn,:));
        end
        
        raw_ax_rf = mean(tmp_ax_rf,1);
        raw_ax_env = mean(tmp_ax_env,1);
        
        ax_res_rf(n,i) = 2*z(find(raw_ax_rf(center_j:end) <= thresh,1));
        ax_res_env(n,i) = 2*z(find(raw_ax_env(center_j:end) <= thresh,1));
        clear tmp_ax_rf tmp_ax_env raw_ax_rf raw_ax_env

    end
end

load('./mat_files/hw5_data_full2.mat');
n = 10;

clear n_r
for n_new = 1:length(scatPerResCell)
    n = n+1;
    for i = 1:length(rndSeed)
        n_r{n_new,i} = find(r{n_new,i}>zposlim(1)+.75*lambda & r{n_new,i}<zposlim(2)-.75*lambda);
        
        n_snr{n_new,i} = find(r{n_new,i}>zposlim(1)+2.5*lambda & r{n_new,i}<zposlim(2)-2.5*lambda);
        
        trim_env = env{n_new,i}(n_r{n_new,i},:);
        trim_rf = rf{n_new,i}(n_r{n_new,i},:);
        env_snr = env{n_new,i}(n_snr{n_new,i},:);
        SNR(n,i) = mean(env_snr(:))/std(env_snr(:));
           
        center_i = size(trim_env,2);
        center_j = size(trim_env,1);
        
%         n_ave_i = 11;
%         n_ave_j = 47;

        corr_i = (center_i-floor(n_ave_i/2)):(center_i+floor(n_ave_i/2));
        corr_j = (center_j-floor(n_ave_j/2)):(center_j+floor(n_ave_j/2));
        
        corr_env{n,i} = normxcorr2(trim_env/max(trim_env(:)), trim_env/max(trim_env(:)));
        corr_rf{n,i} = normxcorr2(trim_rf/max(trim_rf(:)), trim_rf/max(trim_rf(:)));
        
        dz = (1/fs)*(1540/2);
        z = (0:center_j-1).*dz;
        x = linspace(0,max(xposlim),center_i);
        
%         for nn = 1:n_ave_j
%             tmp_lat_rf = abs(hilbert(corr_rf{n,i}(corr_j(nn),:)));
%             tmp_lat_rf = tmp_lat_rf/max(tmp_lat_rf);
%             tmp_lat_env = corr_env{n,i}(corr_j(nn),:);
%             tmp_lat_env = tmp_lat_env/max(tmp_lat_env);
% 
%             raw_lat_rf(nn) = 2*x(find(tmp_lat_rf(center_i:end) <= thresh,1));
%             raw_lat_env(nn) = 2*x(find(tmp_lat_env(center_i:end) <= thresh,1));
%             
% %             keyboard
% 
%             clear tmp_lat_rf tmp_lat_env  
%         end
%         for nn = 1:n_ave_i
%             tmp_ax_rf = abs(hilbert(corr_rf{n,i}(:,corr_i(nn))));
%             tmp_ax_rf = tmp_ax_rf/max(tmp_ax_rf);
%             tmp_ax_env = corr_env{n,i}(:,corr_i(nn));
%             tmp_ax_env = tmp_ax_env/max(tmp_ax_env);
%             
%             raw_ax_rf(nn) = 2*z(find(tmp_ax_rf(center_j:end) <= thresh,1));
%             raw_ax_env(nn) = 2*z(find(tmp_ax_env(center_j:end) <= thresh,1));
%             
%             clear tmp_ax_rf tmp_ax_env 
%         end
%         
%         lat_res_rf(n,i) = mean(raw_lat_rf);
%         lat_res_env(n,i) = mean(raw_lat_env);
%         ax_res_rf(n,i) = mean(raw_ax_rf);
%         ax_res_env(n,i) = mean(raw_ax_env);
%         clear raw_lat_rf raw_lat_env raw_ax_rf raw_ax_env
        
        for nn = 1:n_ave_j
            tmp_lat_rf(nn,:) = abs(hilbert(corr_rf{n,i}(corr_j(nn),:)));
            tmp_lat_rf(nn,:) = tmp_lat_rf(nn,:)/max(tmp_lat_rf(nn,:));
            tmp_lat_env(nn,:) = corr_env{n,i}(corr_j(nn),:);
            tmp_lat_env(nn,:) = tmp_lat_env(nn,:)/max(tmp_lat_env(nn,:));
        end
        
        raw_lat_rf = mean(tmp_lat_rf,1);
        raw_lat_env = mean(tmp_lat_env,1);
        lat_res_rf(n,i) = 2*x(find(raw_lat_rf(center_i:end) <= thresh,1));
        lat_res_env(n,i) = 2*x(find(raw_lat_env(center_i:end) <= thresh,1));
        clear tmp_lat_rf tmp_lat_env raw_lat_rf raw_lat_env
        
        for nn = 1:n_ave_i
            tmp_ax_rf(nn,:) = abs(hilbert(corr_rf{n,i}(:,corr_i(nn))));
            tmp_ax_rf(nn,:) = tmp_ax_rf(nn,:)/max(tmp_ax_rf(nn,:));
            tmp_ax_env(nn,:) = corr_env{n,i}(:,corr_i(nn));
            tmp_ax_env(nn,:) = tmp_ax_env(nn,:)/max(tmp_ax_env(nn,:));
        end
        
        raw_ax_rf = mean(tmp_ax_rf,1);
        raw_ax_env = mean(tmp_ax_env,1);
        
        ax_res_rf(n,i) = 2*z(find(raw_ax_rf(center_j:end) <= thresh,1));
        ax_res_env(n,i) = 2*z(find(raw_ax_env(center_j:end) <= thresh,1));
        clear tmp_ax_rf tmp_ax_env raw_ax_rf raw_ax_env

    end
end

scatPerResCell = [1:15];

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
% print -dpng -r300 ./hw5_speckle_snr.png


figure, hold on;
e1 = errorbar(scatPerResCell,1000*lat_rf_m,1000*lat_rf_s)
e2 = errorbar(scatPerResCell,1000*lat_env_m,1000*lat_env_s,'r')
p1 = plot(scatPerResCell,1000*lat_res.*ones(length(scatPerResCell),1),'k--')
grid on, xlabel('Average Scatterers Per Resolution Cell')
ylabel('Lateral Speckle Size (mm)')
legend([e1 e2 p1],{'RF','Detected','Lateral resolution (\lambdaz/D)'})
hold off
% print -dpng -r300 ./hw5_lat_autocorr_ave.png

figure, hold on;
e1 = errorbar(scatPerResCell,1000*ax_rf_m,1000*ax_rf_s)
e2 = errorbar(scatPerResCell,1000*ax_env_m,1000*ax_env_s,'r')
p1 = plot(scatPerResCell,1000*ax_res.*ones(length(scatPerResCell),1),'k--')
grid on, xlabel('Average Scatterers Per Resolution Cell')
ylabel('Axial Speckle Size (mm)')
legend([e1 e2 p1],{'RF','Detected','Axial resolution (\lambda)'})
hold off
% print -dpng -r300 ./hw5_ax_autocorr_ave.png

