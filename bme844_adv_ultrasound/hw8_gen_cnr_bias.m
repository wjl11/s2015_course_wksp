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

rng(0);

xposlim = [-1 1]*1e-2;
zposlim = [3 5]*1e-2;

maxth = atand(xposlim(2)/foc_z);
th_int = 0.2;
th_scan = (-maxth:th_int:maxth)*pi/180;

scatPerResCell = 5; % minimum for ~ fully developed speckle

scatDensity = scatPerResCell/aResCell;
scatN = round((diff(xposlim)*diff(zposlim))*scatDensity);

biasN = 30;
bias_th = 20+10.*rand(biasN,1);
bias_r = 0.038+0.004.*rand(biasN,1);

bias_x = bias_r.*sind(bias_th);
bias_z = bias_r.*cosd(bias_th);

xpos = [xposlim(1)+(xposlim(2)-xposlim(1)).*rand(scatN,1); bias_x];
zpos = [zposlim(1)+(zposlim(2)-zposlim(1)).*rand(scatN,1); bias_z];
ypos = zeros(scatN+biasN,1);
position = [xpos ypos zpos];

les_amp = ones(scatN,1);
k = find((xpos.^2+(zpos-0.04).^2)<0.005^2);
les_amp(k) = 2;

scat_amp = [1 200];
figure
plot(position(:,1),position(:,3),'x')

for n_amp = 1:length(scat_amp)
    disp(['Simulating off-axis scatterers of amplitude ' num2str(scat_amp(n_amp))])
    amplitude = [les_amp; scat_amp(n_amp).*ones(biasN,1)];

    start=[]; rf =[];
    for nn = 1:length(th_scan)
        xdc_center_focus(tx,[0 0 0]); % Make sure to set center first, then focus
        xdc_focus(tx,0,foc_z.*[sin(th_scan(nn)) 0 cos(th_scan(nn))]);
        xdc_center_focus(rx,[0 0 0]);
        xdc_dynamic_focus(rx,0,th_scan(nn),0);
        [temp,start(nn)]=calc_scat(tx,rx,position,amplitude);
        rf(1:length(temp),nn)=temp;
    end

    min_start=min(start);
    for nn=1:size(rf,2), % shift all RF signals appropriately to generate image
        temp=[zeros(round((start(nn)-min_start)*fs),1);rf(:,nn)];
        temp_rf(1:length(temp),nn)=temp;
    end
    rf=temp_rf; clear temp_rf;

    env=abs(hilbert(rf)); % envelope calculation
    env_db=20*log10(env/max(env(:))); % log scale from 0 to - infinity dB
    r=(min_start+(1:size(env_db,1))/fs) * 1540/2; % determine z axis values
    
    save(['./mat_files/hw8_data' num2str(scat_amp(n_amp)) '_les.mat'],'env','r')
    clear env env_db r
end

save(['./mat_files/hw8_params_les.mat'],'scat_amp','th_scan','position')

xdc_free(tx); xdc_free(rx);
field_end
