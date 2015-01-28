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

for condition = 1:2
    % Create transmit and receive arrays
    if condition == 1
        focus=[0 0 elv_focus]; %focus at 4 cm depth
    elseif condition == 2
        focus=[0 0 elv_focus+0.01];
    end
tx=xdc_linear_array(N_el,el_width,el_height,el_kerf,n_sub_x,n_sub_y,focus);
rx=xdc_linear_array(N_el,el_width,el_height,el_kerf,n_sub_x,n_sub_y,focus);

% Set impulse responses
tc=gauspuls('cutoff',f0,BW); % Note: default BWR -6 dB and TPE -60 dB
imp_resp=gauspuls(-tc:1/fs:tc,f0,BW);
xdc_impulse(tx,imp_resp); 
xdc_impulse(rx,imp_resp);

% positions_pts=[zeros(6,2) (1:6)']*1e-2; %create point targets spaced 1 cm from 1 to 6 cm on z axis
z_array=(15:0.5:80).*1e-3;
x=(-0.8:0.01:0.8).*1e-2;
foc_z=focus(3);

position=[x' zeros(length(x),1) elv_focus*ones(length(x),1)];
rf{condition}=[];start=[];
xdc_center_focus(tx,[0 0 0]); % Make sure to set center first, then focus
xdc_focus(tx,0,[0 0 foc_z]);
xdc_center_focus(rx,[0 0 0]);
xdc_dynamic_focus(rx,0,0,0);
[rf{condition},min_start]=calc_hhp(tx,rx,position);

env{condition}=abs(hilbert(rf{condition})); % envelope calculation
env_db{condition}=20*log10(env{condition}/max(env{condition}(:))); % log scale from 0 to - infinity dB
z{condition}=(min_start+(1:size(env_db{condition},1))/fs) * 1540/2; % determine z axis values
    
clear tmp_i
xdc_free(rx); xdc_free(tx);
clear rx tx
end

% save('hw2_foc_psf_v2.mat');
field_end
%% plot
close all
xlims = [-8 8];
ylims = [38 42];
figure
for i = 1:length(z_array)
    subplot(121)
    imagesc(x*1000,z{1}*1000,rf{1});axis image; colormap gray;
    xlabel('Lateral Position (mm)','Fontsize',14)
    ylabel('Axial Position (mm)','Fontsize',14)
    title('Focused Tx PSF','Fontsize',14)
    ylim([ylims])
    xlim([xlims])

    subplot(122)
    imagesc(x*1000,z{2}*1000,rf{2});axis image; colormap gray;
    xlabel('Lateral Position (mm)','Fontsize',14)
    ylabel('Axial Position (mm)','Fontsize',14)
    title('PW Tx PSF','Fontsize',14)
    ylim([ylims])
    xlim(xlims)
    set(gcf,'position',[100 100 900 300])
end

figure
for i = 1:length(z_array)
    subplot(121)
    imagesc(x*1000,z{1}*1000,env_db{1},[-40 0]);axis image; colormap gray;
    xlabel('Lateral Position (mm)','Fontsize',14)
    ylabel('Axial Position (mm)','Fontsize',14)
    title('Focused Tx PSF','Fontsize',14)
    ylim([ylims])
    xlim([xlims])

    subplot(122)
    imagesc(x*1000,z{2}*1000,env_db{2},[-40 0]);axis image; colormap gray;
    xlabel('Lateral Position (mm)','Fontsize',14)
    ylabel('Axial Position (mm)','Fontsize',14)
    title('PW Tx PSF','Fontsize',14)
    ylim([ylims])
    xlim([xlims])
    set(gcf,'position',[100 100 900 300])
end