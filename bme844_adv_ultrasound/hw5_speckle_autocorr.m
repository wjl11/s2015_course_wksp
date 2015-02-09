close all; clear all; clc;

load('hw5_data.mat');
for n = 1:length(scatN_array)
env_db = 20*log10(env{n}/max(env{n}(:)));
figure;
imagesc(deg_scan,r{n},env{n}); colormap gray;

corr_i = find(th_scan == 0);
corr_j = find(r{n} >= foc_z,1);

[ax_corr_env{n},ax_corr_lag] = xcorr(env{n}(n_r{n},corr_i),'coeff');
ax_corr_rf{n} = xcorr(rf{n}(n_r{n},corr_i),'coeff');

ax_lag{n} = ax_corr_lag.*(1/fs)*1540/2;

[lat_corr_env{n},lat_corr_lag] = xcorr(env{n}(corr_j,:),'coeff');
lat_corr_rf{n} = xcorr(rf{n}(corr_j,:),'coeff');

lat_lag{n} = lat_corr_lag.*th_int;

corr_env{n} = xcorr2(env{n}/max(env{n}(:)));
corr_rf{n} = xcorr2(rf{n}/max(rf{n}(:)));

% 
% figure;
% imagesc(deg_scan(n_deg),r{n}(n_r{n}),env_db(n_r{n},n_deg),[-40 0]); colormap gray;

figure;
hold on; p1 = plot(lat_lag{n},lat_corr_env{n},'r'); p2 = plot(lat_lag{n}, lat_corr_rf{n},'b'); hold on;
legend([p1,p2],{'Detected','RF'}),ylabel('Lateral autocorrelation'); xlabel('Lag (cm)')

figure; hold on; plot(ax_lag{n},ax_corr_env{n},'r'); plot(ax_lag{n}, ax_corr_rf{n},'b'); hold on;

end