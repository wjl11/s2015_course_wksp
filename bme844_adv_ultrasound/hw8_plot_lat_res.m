close all; clear all; clc;
load('./mat_files/hw8_params.mat')
scat_amp = scat_amp;
idx = 1;
figure
for i = scat_amp
    load(['./mat_files/hw8_data' num2str(i)])
    env_db=20*log10(env/max(env(:))); % log scale from 0 to - infinity dB
    
    imagesc(1000*x_scan,1000*r,env_db,[-40 0]); colormap gray
    xlabel('Lateral Position (mm)'), ylabel('Depth (mm)');
    title(['PSF at 0^o (off-axis scatterering amplitude ' num2str(i) ')']);
    
    c6tmp = contourc(x_scan,r,env_db,[-6 -6]);
    c12tmp = contourc(x_scan,r,env_db,[-12 -12]);

    tmp_i=find(c6tmp(1,:)~=-6);
    c6db=c6tmp(:,tmp_i);
    tmp_i = find(c12tmp(1,:)~=-12);
    c12db=c12tmp(:,tmp_i);
    
    fw6db(idx)=max(c6db(1,:))-min(c6db(1,:));
    fw12db(idx)=max(c12db(1,:))-min(c12db(1,:));
    
    axfw6db(idx)=max(c6db(2,:))-min(c6db(2,:));
    axfw12db(idx)=max(c12db(2,:))-min(c12db(2,:));
    idx = idx+1;
    clear env env_db r c6db c12db
    pause
end

figure
plot(scat_amp,fw6db,'x')
xlabel('Scattering Amplitude'); ylabel('Lateral Resolution (6 dB)')
xlim([0 200])

% figure 
% plot(scat_amp,axfw6db,'x')