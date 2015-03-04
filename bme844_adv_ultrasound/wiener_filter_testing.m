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
for n = length(scatPerResCell)
    for i = length(rndSeed)
        n_r{n,i} = find(r{n,i}>zposlim(1)+2*lambda & r{n,i}<zposlim(2)-2*lambda);
        
        n_snr{n,i} = find(r{n,i}>zposlim(1)+2.5*lambda & r{n,i}<zposlim(2)-2.5*lambda);
        
        trim_env = env{n,i}(n_r{n,i},:);
        trim_rf = rf{n,i}(n_r{n,i},:);
        env_snr = env{n,i}(n_snr{n,i},:);
        SNR(n,i) = mean(env_snr(:))/std(env_snr(:));
           
        center_i = size(trim_env,2);
        center_j = size(trim_env,1);
        
        corr_i = (center_i-floor(n_ave_i/2)):(center_i+floor(n_ave_i/2));
        corr_j = (center_j-floor(n_ave_j/2)):(center_j+floor(n_ave_j/2));
        
        corr_env = normxcorr2(trim_env/max(trim_env(:)), trim_env/max(trim_env(:)));
        corr_rf = normxcorr2(trim_rf/max(trim_rf(:)), trim_rf/max(trim_rf(:)));
        
        trim_env = trim_env/max(trim_env(:));
        trim_rf = trim_rf/max(trim_rf(:));
        
        rf_PSF = corr_rf(150:300,60:140);
        
        wnr1 = deconvwnr(trim_rf, rf_PSF, 0);
        
        wnr1env = abs(hilbert(wnr1));
        trim_env = abs(hilbert(trim_rf));
        
        figure
        subplot(121)
        imagesc(20*log10(wnr1env/max(wnr1env(:))) ,[-40 0]); colormap gray;
        subplot(122)
        imagesc(20*log10(trim_env/max(trim_env(:))) ,[-40 0]); colormap gray;
%         
%         dz = (1/fs)*(1540/2);
%         z = (0:center_j-1).*dz;
%         x = linspace(0,max(xposlim),center_i);
    end
end

