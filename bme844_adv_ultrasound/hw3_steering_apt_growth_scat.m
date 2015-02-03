close all; clear all; clc;
addpath('../../field_ii')
field_init(-1)
% Transducer parameters
f0=5e6; BW=0.7; 
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

maxth = 10;
th_scan = [-10:0.1:maxth]*pi/180;
th_target = [0:5:maxth-5]*pi/180;
foc_z = focus(3);

N_el_tot = 128;
constFnum = foc_z/((N_el_tot*el_width)*cosd(maxth))

for i = 1:length(th_target)

%     if th_target(i) == maxth
%         N_el = N_el_tot;
%     else
%         tmp = foc_z/(el_width*constFnum*cos(th_target(i)));
%         N_el = ceil(tmp);
%         clear tmp
%         N_el = N_el_tot;
%     end
%     Fnum(i) = foc_z/(N_el*el_width*cos(th_target(i)));
    
    tx=xdc_linear_array(N_el_tot,el_width,el_height,el_kerf,n_sub_x,n_sub_y,focus);
    rx=xdc_linear_array(N_el_tot,el_width,el_height,el_kerf,n_sub_x,n_sub_y,focus);

    % Set impulse responses
    tc=gauspuls('cutoff',f0,BW); % Note: default BWR -6 dB and TPE -60 dB
    imp_resp=gauspuls(-tc:1/fs:tc,f0,BW);
    xdc_impulse(tx,imp_resp); 
    xdc_impulse(rx,imp_resp);

    start=[];
    position = foc_z.*[sin(th_target(i)) 0 cos(th_target(i))];
    amplitude = 100;
    for nn = 1:length(th_scan)
        
        if th_scan(nn) == deg2rad(maxth)
            apod = ones(1,N_el_tot);
        else
            tmp1 = foc_z/(constFnum*el_width*cos(th_scan(nn)));
            tmp2 = ones(1,round(tmp1));
            pad = floor((N_el_tot-length(tmp2))/2);
           	fill = N_el_tot-2*pad-length(tmp2);
            apod = [padarray(tmp2,[0,pad]) zeros(1,fill)];
        end
        el(nn) = length(find(apod == 1));
        Fnum(nn) = foc_z/(el(nn)*el_width*cos(th_scan(nn)));
        
        xdc_apodization(tx,0,apod);
        xdc_apodization(rx,0,apod);
        xdc_center_focus(tx,[0 0 0]); % Make sure to set center first, then focus
        xdc_focus(tx,0,foc_z.*[sin(th_scan(nn)) 0 cos(th_scan(nn))]);
        xdc_center_focus(rx,[0 0 0]);
%         xdc_focus(rx,0,foc_z.*[sin(th(nn)) 0 cos(th(nn))]);
        xdc_dynamic_focus(rx,0,th_scan(nn),0);
        [temp,start(nn)]=calc_scat(tx,rx,position,amplitude);
        rf{i}(1:length(temp),nn)=temp;
    end
    
    min_start=min(start);
    for nn=1:size(rf{i},2), % shift all RF signals appropriately to generate image
        temp=[zeros(round((start(nn)-min_start)*fs),1);rf{i}(:,nn)];
        temp_rf(1:length(temp),nn)=temp;
    end
    rf{i}=temp_rf; clear temp_rf;
    
    env{i}=abs(hilbert(rf{i})); % envelope calculation
    env_db=20*log10(env{i}/max(env{i}(:))); % log scale from 0 to - infinity dB
    r{i}=(min_start+(1:size(env_db,1))/fs) * 1540/2; % determine z axis values
    
%     [out,ax,l,outside] = scan_convert(env_db{i},'sector',0,maxth,0,2,fs);
    
    idx = find(r{i}*100>3.8 & r{i}*100<4.2);
    env_trim{i} = env{i}(idx,:);
    rf_trim{i} = rf{i}(idx,:);
    
%     figure
% %     imagesc(l,ax,out,[-40 0]); axis image; colormap gray;
%     imagesc(th_scan,r{i}*100,env_trim{i},[-40 0]); axis image; colormap jet;
% 
%     figure
%     imagesc(th_scan,r{i}(idx)*100,rf_trim{i}); axis image; colormap jet;
%     
    
    C_env{i} = xcorr2(env_trim{1},env_trim{i});
    C_rf{i} = xcorr2(rf_trim{1},rf_trim{i});
    
    figure
    imagesc(th_scan,r{i}*100,C_env{i}); axis image; colormap jet; colorbar;
    
    figure
    imagesc(th_scan,r{i}*100,C_rf{i}); axis image; colormap jet; colorbar;
    
    xdc_free(tx); xdc_free(rx);
end

field_end


