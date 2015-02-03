close all; clear all; clc;
addpath('../../field_ii')
field_init(-1)
% Transducer parameters
f0=5e6; BW=0.7; N_el=128; %5 MHz, 128 element linear array w/ lamda pitch, 70% BW, 6.5 F number
elv_focus=0.04; elv_fnum=6.5; kerf_fraction=0.05;

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

tx=xdc_linear_array(N_el,el_width,el_height,el_kerf,n_sub_x,n_sub_y,focus);
rx=xdc_linear_array(N_el,el_width,el_height,el_kerf,n_sub_x,n_sub_y,focus);

% Set impulse responses
tc=gauspuls('cutoff',f0,BW); % Note: default BWR -6 dB and TPE -60 dB
imp_resp=gauspuls(-tc:1/fs:tc,f0,BW);
xdc_impulse(tx,imp_resp); 
xdc_impulse(rx,imp_resp);

th = [0:10:40]*pi/180;
foc_z = focus(3);

% z_array = foc_z.*cosd(ang);
% x_array = foc_z.*sind(ang);
% xpos=(0:0.01:2).*1e-2;
for i = 1:length(th)
    if th(i) > 0
        zpos = (xpos/tan(th(i)));
    else
        zpos = foc_z*ones(1,length(xpos));
    end
    position = [xpos' zeros(length(xpos),1) zpos'];

    xdc_center_focus(tx,[0 0 0]); % Make sure to set center first, then focus
    xdc_focus(tx,0,foc_z.*[sin(th(i)) 0 cos(th(i))]);
    xdc_center_focus(rx,[0 0 0]);
    xdc_dynamic_focus(rx,0,th(i),0);
    
%     [rf{i},min_start]=calc_hhp(tx,rx,position);
    
    env=abs(hilbert(rf{i})); % envelope calculation
    env_db{i}=20*log10(env/max(env(:))); % log scale from 0 to - infinity dB
    z{i}=(min_start+(1:size(env_db{i},1))/fs) * 1540/2; % determine z axis values
    
    figure
    imagesc(xpos*100,z{i}*100,env_db{i},[-40 0]);axis image; colormap gray;
    
%     xdc_free(tx); xdc_free(rx);
end
field_end


