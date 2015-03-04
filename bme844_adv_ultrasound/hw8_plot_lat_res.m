close all; clear all; clc;
load('./mat_files/hw8_params.mat')
scat_amp = scat_amp;
idx = 1;
figure

writerObj1 = VideoWriter('hw8_offaxis_scattering.avi')
writerObj1.FrameRate = 1;
open(writerObj1);

for i = scat_amp
    load(['./mat_files/hw8_data' num2str(i)])
    env_db=20*log10(env/max(env(:))); % log scale from 0 to - infinity dB
    trimi = find(1000*x_scan<3 & 1000*x_scan>-3);
    trimj = find(1000*r<42 & 1000*r>39);
    
    imagesc(1000*x_scan(trimi),1000*r(trimj),env_db(trimj,trimi),[-40 0]); colormap gray
    
    xlabel('Lateral Position (mm)'), ylabel('Depth (mm)');
    title(sprintf('PSF at 0^o (off-axis scattering amplitude %.1f dB)',20.*log10(i)));
    frame = getframe(gcf);
    writeVideo(writerObj1,frame);
    
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
end
close(writerObj1)

figure
plot(scat_amp,1000.*fw6db,'x')
xlabel('Scattering Amplitude'); ylabel('6 dB Lateral Resolution (mm)')
xlim([0 200])
print -dpng -r300 ./hw8_lat_res_psf.png
