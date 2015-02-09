close all; clear all; clc;
addpath('../../field_ii')
field_init(-1)

% Transducer parameters
f0=5e6; BW=0.75; 
% N_el=128; %5 MHz, 128 element linear array w/ lamda pitch, 70% BW, 6.5 F number
elv_focus=0.04; elv_fnum=6.5; kerf_fraction=0;
% NOTE: calculations for fnumber do not account for kerf_fraction (keep value at 0)

% General parameters
c=1540; fs=100e6; % speed of sound c and ultrasound sampling f (default fs)

% Derived parameters
lambda=c/f0;
pitch=lambda/2;
el_height=elv_focus/elv_fnum; % element height (size in y direction)
el_width=(1-kerf_fraction)*pitch; % element width (x direction)
el_kerf=kerf_fraction*pitch;
n_sub_x=ceil(el_width/(lambda/4));
n_sub_y=ceil(el_height/(lambda/4));

focus=[0 0 elv_focus]; %focus at 4 cm depth

maxth = 5;
th_int = 0.1;
th_scan = [-maxth:th_int:maxth]*pi/180;
foc_z = focus(3);

N_el_tot = 128;
tx=xdc_linear_array(N_el_tot,el_width,el_height,el_kerf,n_sub_x,n_sub_y,focus);
rx=xdc_linear_array(N_el_tot,el_width,el_height,el_kerf,n_sub_x,n_sub_y,focus);

D = N_el_tot*(el_width+el_kerf)-el_kerf;

% Set impulse responses
tc=gauspuls('cutoff',f0,BW); % Note: default BWR -6 dB and TPE -60 dB
imp_resp=gauspuls(-tc:1/fs:tc,f0,BW);
xdc_impulse(tx,imp_resp); 
xdc_impulse(rx,imp_resp);

lat_res = lambda*foc_z/D;
ax_res = lambda;
aResCell = lat_res*ax_res;

ncell = 10;
xposlim = [-ncell*lat_res ncell*lat_res];
zposlim = [foc_z-ncell*ax_res foc_z+ncell*ax_res];

scatN_array = 4.*[100 300 500 700];

for n = 1:length(scatN_array)
scatN = scatN_array(n);
xpos = xposlim(1)+(xposlim(2)-xposlim(1)).*rand(scatN,1);
zpos = zposlim(1)+(zposlim(2)-zposlim(1)).*rand(scatN,1);
ypos = zeros(scatN,1);
scatDensity(n) = scatN/(diff(xposlim)*diff(xposlim));
scatPerResCell(n) = scatDensity(n)*aResCell
position = [xpos ypos zpos];
amplitude = ones(scatN,1);

start=[]; 
for nn = 1:length(th_scan)
    xdc_center_focus(tx,[0 0 0]); % Make sure to set center first, then focus
    xdc_focus(tx,0,foc_z.*[sin(th_scan(nn)) 0 cos(th_scan(nn))]);
    xdc_center_focus(rx,[0 0 0]);
    xdc_dynamic_focus(rx,0,th_scan(nn),0);
    [temp,start(nn)]=calc_scat(tx,rx,position,amplitude);
    rf{n}(1:length(temp),nn)=temp;
    disp(['Scan line: ' num2str(nn) '/' num2str(length(th_scan))])
end

min_start=min(start);
for nn=1:size(rf{n},2), % shift all RF signals appropriately to generate image
    temp=[zeros(round((start(nn)-min_start)*fs),1);rf{n}(:,nn)];
    temp_rf(1:length(temp),nn)=temp;
end
rf{n}=temp_rf; clear temp_rf;

env{n}=abs(hilbert(rf{n})); % envelope calculation
env_db=20*log10(env{n}/max(env{n}(:))); % log scale from 0 to - infinity dB
r{n}=(min_start+(1:size(env_db,1))/fs) * 1540/2; % determine z axis values

deg_scan = rad2deg(th_scan);
n_deg = find(deg_scan>-4 & deg_scan<4);
n_r{n} = find(r{n}>0.0395 & r{n}<0.042);

% corr_i = find(th_scan == 0);
% corr_j = find(r{n} >= foc_z,1);
% 
% [ax_corr_env{n},ax_corr_lag] = autocorr(env{n}(n_r{n},corr_i),round(size(env{n},1)/2)-1);
% ax_corr_rf{n} = autocorr(rf{n}(n_r{n},corr_i),round(size(env{n},1)/2)-1);
% 
% ax_lag{n} = ax_corr_lag.*(1/fs)*1540/2;
% 
% [lat_corr_env{n},lat_corr_lag] = autocorr(env{n}(corr_j,n_deg),round(size(env{n},2)/2)-1);
% lat_corr_rf{n} = autocorr(rf{n}(corr_j,n_deg),round(size(env{n},2)/2)-1);
% 
% lat_lag{n} = lat_corr_lag.*th_int;
% 
% corr_env{n} = xcorr2(env{n}(n_r{n},n_deg));
% corr_rf{n} = xcorr2(rf{n}(n_r{n},n_deg));
% 
figure;
imagesc(deg_scan,r{n},env_db,[-40 0]); colormap gray;

% figure;
% hold on; plot(lat_lag{n},lat_corr_env{n},'r'); plot(lat_lag{n}, lat_corr_rf{n},'b'); hold on;
% figure; hold on; plot(ax_lag{n},ax_corr_env{n},'r'); plot(ax_lag{n}, ax_corr_rf{n},'b'); hold on;

end

save('hw5_data.mat')
xdc_free(tx); xdc_free(rx);
field_end