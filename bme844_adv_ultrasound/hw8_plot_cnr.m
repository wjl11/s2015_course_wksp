clear all; close all; clc;

load('./mat_files/hw8_params_les.mat');

for n = 1:length(scat_amp)
        filename = ['./mat_files/hw8_data' num2str(scat_amp(n)) '_les.mat'];
        disp(['Loading ' filename]);
        load(filename);

        env_trim = env(1100:3200,10:130);
        env_mat(:,:,n) = env_trim./mean(env_trim(:));
end

th_les = [60:80];
r_les = [600:1200];

th_bg = [100:120];
r_bg = [1000:1600];

for n = 1:length(scat_amp)

    im_les = env_mat(r_les,th_les,n);
    im_bg = env_mat(r_bg,th_bg,n);
    
    tmp = im_les(:);
    u_les = mean(tmp(~isnan(tmp)));
    v_les = var(tmp(~isnan(tmp)));
    clear tmp

    tmp = im_bg(:);
    u_bg = mean(tmp(~isnan(tmp)));
    v_bg = var(tmp(~isnan(tmp)));
    clear tmp

    CNR(n) = abs(u_bg-u_les)/sqrt(v_bg+v_les);
    
    clear im_les im_bg u_bg v_bg u_les v_les
end

figure
plot(scat_amp(1:end-1),CNR(1:end-1),'x')
xlabel('Scattering Amplitude'),ylabel('CNR')
print -dpng -r300 ./hw8_CNR.png

