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
        focus=[0 0 10];
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

for i = 1:length(z_array)
    position=[x' zeros(length(x),1) z_array(i)*ones(length(x),1)];
    rf=[];start=[];
    xdc_center_focus(tx,[0 0 0]); % Make sure to set center first, then focus
    xdc_focus(tx,0,[0 0 foc_z]);
    xdc_center_focus(rx,[0 0 0]);
    xdc_dynamic_focus(rx,0,0,0);
    [rf,min_start]=calc_hhp(tx,rx,position);
    
    env=abs(hilbert(rf)); % envelope calculation
    env_db{condition,i}=20*log10(env/max(env(:))); % log scale from 0 to - infinity dB
    z{condition,i}=(min_start+(1:size(env_db{condition,i},1))/fs) * 1540/2; % determine z axis values

    c6tmp{condition,i} = contourc(x*100,z{condition,i}*100,env_db{condition,i},[-6 -6]);
    c12tmp{condition,i} = contourc(x*100,z{condition,i}*100,env_db{condition,i},[-12 -12]);

    tmp_i=find(c6tmp{condition,i}(1,:)~=-6);
    c6db{condition,i}=c6tmp{condition,i}(:,tmp_i);
    tmp_i = find(c12tmp{condition,i}(1,:)~=-12);
    c12db{condition,i}=c12tmp{condition,i}(:,tmp_i);
    
    fw6db(i,condition)=max(c6db{condition,i}(1,:))-min(c6db{condition,i}(1,:));
    fw12db(i,condition)=max(c12db{condition,i}(1,:))-min(c12db{condition,i}(1,:));
    disp(sprintf('%.2f cm depth processed.',z_array(i)*100))
end
clear env tmp_i rf
xdc_free(rx); xdc_free(tx);
clear rx tx
end

save('sim_data_hhp.mat');
field_end
%% plot
close all
figure
for i = 1:length(z_array)
    subplot(121)
    imagesc(x*100,z{1,i}*100,env_db{1,i},[-40 0]);axis image; colormap gray;
%     t1 = title(sprintf('Focused Tx (z = %.2f cm)',100*z_array(i)),'fontsize',16,'Color','w')
    uicontrol('Style', 'text',...
       'String', 'Focused Tx',... %replace something with the text you want
       'Units','normalized',...
       'ForegroundColor','w',...  
       'BackgroundColor','k',...
       'Fontsize',18,...
       'Position', [0.2 0.75 0.2 0.07]); 
   
%     set(t1,'position',[250 300 0])
    axis off;
    
    uicontrol('Style', 'text',...
       'String',  sprintf('z = %.2f cm',100*z_array(i)),... %replace something with the text you want
       'Units','normalized',...
       'ForegroundColor','w',...  
       'BackgroundColor','k',...
       'Fontsize',24,...
       'Position', [0.4 0.25 0.2 0.07]); 
    
    subplot(122)
    imagesc(x*100,z{2,i}*100,env_db{2,i},[-40 0]);axis image; colormap gray;
%     t2 = title(,'fontsize',16,'Color','w')
    uicontrol('Style', 'text',...
       'String', 'Unfocused Tx',... %replace something with the text you want
       'Units','normalized',...
       'ForegroundColor','w',...  
       'BackgroundColor','k',...
       'Fontsize',18,...
       'Position', [0.65 0.75 0.2 0.07]); 
    set(gcf,'color','k');
%     set(t2,'position',[750 300 0])
    axis off
    set(gcf,'position',[100 100 1000 400])
    pause(0.1)
end