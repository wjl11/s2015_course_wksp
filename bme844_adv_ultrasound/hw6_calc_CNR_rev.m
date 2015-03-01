clear all; close all; clc;

load('hw6_params_test.mat');

df = [0:0.2:1.2];
for n = 1:length(df)
    idx = 1;
    if df(n) == 0 
        filename = ['hw6_data8.mat'];
        disp(['Loading ' filename]);
        load(filename);
        
        env_trim = env(1100:3200,10:130);
        env_mat(:,:,idx) = env_trim./max(env_trim(:));
        env_cont = env_mat(:,:,idx);
        imagesc(env_cont)
        clear env_trim
    else
        for f = 2:df(n):df(n)*floor(10/df(n))
            filename = ['hw6_data' num2str(f) '.mat'];
            disp(['Loading ' filename]);
            load(filename);

            env_trim = env(1100:3200,10:130);
            env_mat(:,:,idx) = env_trim./max(env_trim(:));
            idx = idx+1;
            clear env_trim
        end
    end
    mean_env{n} = mean(env_mat,3);
    clear env_mat idx
end

th_les = [60:80];
r_les = [600:1200];

th_bg = [100:120];
r_bg = [1000:1600];

for n = 1:length(df)

    im_les = mean_env{n}(r_les,th_les);
    im_bg = mean_env{n}(r_bg,th_bg);
    
    tmp = im_les(:);
    u_les = mean(tmp(~isnan(tmp)));
    v_les = var(tmp(~isnan(tmp)));
    clear tmp

    tmp = im_bg(:);
    u_bg = mean(tmp(~isnan(tmp)));
    v_bg = var(tmp(~isnan(tmp)));
    clear tmp

    CNR(n) = abs(u_bg-u_les)/sqrt(v_bg+v_les);
    
%     figure
%     subplot(211)
%     imagesc(im_les); colormap gray;
%     subplot(212)
%     imagesc(im_bg); colormap gray;
    
    clear im_les im_bg u_bg v_bg u_les v_les
end


env_maxCNR = 20*log10(mean_env{3}/max(mean_env{3}(:)));
env_conCNR = 20*log10(env_cont/max(env_cont(:)));

figure
imagesc(env_maxCNR,[-40 0]); colormap gray; axis square
xlabel('Angle (degrees)'),ylabel('Depth (cm)')
% print -dpng -r300 ./hw6_compound.png

figure 
imagesc(env_conCNR,[-40 0]); colormap gray; axis square
xlabel('Angle (degrees)'),ylabel('Depth (cm)')
% print -dpng -r300 ./hw6_control.png

CNR_ratio = CNR(2:end)./CNR(1);

figure
plot(df(2:end),CNR_ratio)
xlabel('\Deltaf (MHz)'),ylabel('Contrast Ratio (CNR_a_v_g/CNR_o_r_i_g)')
% print -dpng -r300 ./hw6_CNR.png

