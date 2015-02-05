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

maxth = 60;
th_scan = [-5:0.1:maxth]*pi/180;
th_target = [0:10:maxth-10]*pi/180;

% maxth = 70;
% th_scan = [-5:0.1:maxth]*pi/180;
% th_target = [0:5:maxth-10]*pi/180;
foc_z = focus(3);

N_el_tot = 128;
constFnum = foc_z/((N_el_tot*el_width)*cosd(maxth));

for n = 1:2
    tx=xdc_linear_array(N_el_tot,el_width,el_height,el_kerf,n_sub_x,n_sub_y,focus);
    rx=xdc_linear_array(N_el_tot,el_width,el_height,el_kerf,n_sub_x,n_sub_y,focus);

    % Set impulse responses
    tc=gauspuls('cutoff',f0,BW); % Note: default BWR -6 dB and TPE -60 dB
    imp_resp=gauspuls(-tc:1/fs:tc,f0,BW);
    xdc_impulse(tx,imp_resp); 
    xdc_impulse(rx,imp_resp);
    for i = 1:length(th_target)
        start=[];
        position = foc_z.*[sin(th_target(i)) 0 cos(th_target(i))];
        amplitude = 100;
        for nn = 1:length(th_scan)
            if n == 1
                if th_scan(nn) == deg2rad(maxth)
                    apod = ones(1,N_el_tot);
                else
                    tmp1 = foc_z/(constFnum*el_width*cos(th_scan(nn)));
                    tmp2 = ones(1,round(tmp1));
                    pad = floor((N_el_tot-length(tmp2))/2);
                    fill = N_el_tot-2*pad-length(tmp2);
                    apod = [padarray(tmp2,[0,pad]) zeros(1,fill)];
                end
                clear tmp1 tmp2
                el(nn) = length(find(apod == 1));
                Fnum(nn) = foc_z/(el(nn)*el_width*cos(th_scan(nn)));

                xdc_apodization(tx,0,apod);
                xdc_apodization(rx,0,apod);
            else
                xdc_apodization(tx,0,ones(1,N_el_tot));
                xdc_apodization(rx,0,ones(1,N_el_tot));
            end
            
            xdc_center_focus(tx,[0 0 0]); % Make sure to set center first, then focus
            xdc_focus(tx,0,foc_z.*[sin(th_scan(nn)) 0 cos(th_scan(nn))]);
            xdc_center_focus(rx,[0 0 0]);
    %         xdc_focus(rx,0,foc_z.*[sin(th(nn)) 0 cos(th(nn))]);
            xdc_dynamic_focus(rx,0,th_scan(nn),0);
            [temp,start(nn)]=calc_scat(tx,rx,position,amplitude);
            rf{n,i}(1:length(temp),nn)=temp;
        end

        min_start=min(start);
        for nn=1:size(rf{n,i},2), % shift all RF signals appropriately to generate image
            temp=[zeros(round((start(nn)-min_start)*fs),1);rf{n,i}(:,nn)];
            temp_rf(1:length(temp),nn)=temp;
        end
        rf{n,i}=temp_rf; clear temp_rf;

        env{n,i}=abs(hilbert(rf{n,i})); % envelope calculation
        env_db=20*log10(env{n,i}/max(env{n,i}(:))); % log scale from 0 to - infinity dB
        r{n,i}=(min_start+(1:size(env_db,1))/fs) * 1540/2; % determine z axis values

        env{n,i} = env{n,i}./max(env{n,i}(:));
        rf{n,i} = rf{n,i}./max(rf{n,i}(:));

    %     [out,ax,l,outside] = scan_convert(env_db{i},'sector',0,maxth,0,2,fs);

        idx = find(r{n,i}*100>3.8 & r{n,i}*100<4.2);
        env_trim{n,i} = env{n,i}(idx,:);
        rf_trim{n,i} = rf{n,i}(idx,:);
        
        tmp1 = diff(r{n,i});
        tmp2 = diff(th_scan);
        T_y(n,i) = 1000*tmp1(1);
        T_x(n,i) = tmp2(1);
        clear tmp1 tmp2
        k_y{n,i} = (1/(2*T_y(n,i)))*linspace(-1,1,length(idx));
        k_x{n,i} = (1/(2*T_x(n,i)))*linspace(-1,1,length(th_scan));

        env_trim_k{n,i} = fftshift(fft2(env_trim{n,i}));
        rf_trim_k{n,i} = fftshift(fft2(rf_trim{n,i}));

        % calculate correlation via xcorr in spatial domain with kernel
        if th_target(i) == 0
            k_i = find(r{n,i}*100>4 & r{n,i}*100<4.15);
            k_j = find(th_scan>-0.05 & th_scan<0.05);
            env_kernel = env{n,i}(k_i,k_j);
            rf_kernel = rf{n,i}(k_i,k_j);
        end

        C_env{n,i} = normxcorr2(env_kernel,env_trim{n,i});
        C_rf{n,i} = normxcorr2(rf_kernel,rf_trim{n,i});

        Cmax_env(n,i) = max(C_env{n,i}(:));
        Cmax_rf(n,i) = max(C_rf{n,i}(:));

        % calculate correlation via normalized k-space product
        kC_env(n,i) = sum(sum(abs(env_trim_k{n,1}).*abs(env_trim_k{n,i})))/...
            sqrt(sum(sum(abs(env_trim_k{n,1}).^2))*sum(sum(abs(env_trim_k{n,i}).^2)));
        kC_rf(n,i) = sum(sum(abs(rf_trim_k{n,1}).*abs(rf_trim_k{n,i})))/...
            sqrt(sum(sum(abs(rf_trim_k{n,1}).^2))*sum(sum(abs(rf_trim_k{n,i}).^2)));
        
        disp(['Image at angle ' num2str(rad2deg(th_target(i))) ' simulated.'])
    end
    xdc_free(tx); xdc_free(rx);
end

figure 
subplot(211), plot(rad2deg(th_scan),el)
grid on
xlabel('Angle (deg)','fontsize',12), ylabel('No. Active Elements','fontsize',12)
subplot(212), plot(rad2deg(th_scan),Fnum)
grid on
xlabel('Angle (deg)','fontsize',12), ylabel('F-Number','fontsize',12)

print -dpng -r300 ./hw4_aperture.png

figure, hold on;
p1 = plot(rad2deg(th_target), Cmax_env(1,:),'bx-','linewidth',2)
p2 = plot(rad2deg(th_target), Cmax_env(2,:),'rx-','linewidth',2)
grid on,ylim([0.90 1])
hold off
legend([p1 p2],{'Aperture growth','Fixed aperture'},'fontsize',12,'location','southwest')
xlabel('Angle (deg)','fontsize',12)
ylabel('Correlation','fontsize',12)
title('PSF Decorrelation (Envelope Detected)','fontsize',12)

print -dpng -r300 ./hw4_psf_corr_env.png

figure, hold on;
p3 = plot(rad2deg(th_target), Cmax_rf(1,:),'bx-','linewidth',2)
p4 = plot(rad2deg(th_target), Cmax_rf(2,:),'rx-','linewidth',2)
grid on, ylim([0.90 1])
hold off
legend([p3 p4],{'Aperture growth','Fixed aperture'},'fontsize',12,'location','southwest')
xlabel('Angle (deg)','fontsize',12)
ylabel('Correlation')
title('PSF Decorrelation (RF)','fontsize',12)

print -dpng -r300 ./hw4_psf_corr_rf.png

figure, hold on;
p1 = plot(rad2deg(th_target), kC_env(1,:),'bx-','linewidth',2)
p2 = plot(rad2deg(th_target), kC_env(2,:),'rx-','linewidth',2)
grid on,ylim([0.90 1])
hold off
legend([p1 p2],{'Aperture growth','Fixed aperture'},'fontsize',12,'location','southwest')
xlabel('Angle (deg)','fontsize',12)
ylabel('Correlation','fontsize',12)
title('PSF Decorrelation (Envelope Detected)','fontsize',12)

print -dpng -r300 ./hw4_psf_corr_env_kspace.png

figure, hold on;
p3 = plot(rad2deg(th_target), kC_rf(1,:),'bx-','linewidth',2)
p4 = plot(rad2deg(th_target), kC_rf(2,:),'rx-','linewidth',2)
grid on, ylim([0.90 1])
hold off
legend([p3 p4],{'Aperture growth','Fixed aperture'},'fontsize',12,'location','southwest')
xlabel('Angle (deg)','fontsize',12)
ylabel('Correlation')
title('PSF Decorrelation (RF)','fontsize',12)

print -dpng -r300 ./hw4_psf_corr_rf_kspace.png


ylims = [-10 10];
xlims = [-100 100];

figure
subplot(131)
imagesc(k_x{1,1},k_y{1,1},abs(rf_trim_k{1,1})); colormap jet; 
ylim(ylims), xlim(xlims)
xlabel('k\theta (1/rads)','fontsize',12), ylabel('kr (1/mm)','fontsize',12)
title(['K-space ' num2str(rad2deg(th_target(1))) ' deg (Aperture Growth)'],'fontsize',12)
subplot(132)
imagesc(k_x{1,2},k_y{1,2},abs(rf_trim_k{1,2})); colormap jet; 
ylim(ylims), xlim(xlims)
xlabel('k\theta (1/rads)','fontsize',12), ylabel('kr (1/mm)','fontsize',12)
title(['K-space ' num2str(rad2deg(th_target(2))) ' deg (Aperture Growth)'],'fontsize',12)
subplot(133)
imagesc(k_x{1,end},k_y{1,end},abs(rf_trim_k{1,end})); colormap jet; 
ylim(ylims), xlim(xlims)
xlabel('k\theta (1/rads)','fontsize',12), ylabel('kr (1/mm)','fontsize',12)
title(['K-space ' num2str(rad2deg(th_target(end))) ' deg (Aperture Growth)'],'fontsize',12)

figure
subplot(131)
imagesc(k_x{2,1},k_y{2,1},abs(rf_trim_k{2,1})); colormap jet;
ylim(ylims), xlim(xlims)
xlabel('k\theta (1/rads)','fontsize',12), ylabel('kr (1/mm)','fontsize',12)
title(['K-space ' num2str(rad2deg(th_target(1))) ' deg (Fixed Aperture)'],'fontsize',12)
subplot(132)
imagesc(k_x{2,2},k_y{2,2},abs(rf_trim_k{2,2})); colormap jet; 
ylim(ylims), xlim(xlims)
xlabel('k\theta (1/rads)','fontsize',12), ylabel('kr (1/mm)','fontsize',12)
title(['K-space ' num2str(rad2deg(th_target(2))) ' deg (Fixed Aperture)'],'fontsize',12)
subplot(133)
imagesc(k_x{2,end},k_y{2,end},abs(rf_trim_k{2,end})); colormap jet; 
ylim(ylims), xlim(xlims)
xlabel('k\theta (1/rads)','fontsize',12), ylabel('kr (1/mm)','fontsize',12)
title(['K-space ' num2str(rad2deg(th_target(end))) ' deg (Fixed Aperture)'],'fontsize',12)

field_end




