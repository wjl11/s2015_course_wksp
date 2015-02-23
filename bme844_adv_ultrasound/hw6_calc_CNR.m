clear all; close all; clc;

load('hw6_params_test.mat');


df = [0:0.2:1.2];
for n = 1:length(df)
    idx = 1;
    if df(n) == 0 
        filename = ['hw6_data5.mat'];
        disp(['Loading ' filename]);
        load(filename);
        env_mat(:,:,idx) = env;
        env_cont = env;
    else
        for f = 2:df(n):df(n)*floor(10/df(n))
            filename = ['hw6_data' num2str(f) '.mat'];
            disp(['Loading ' filename]);
            load(filename);
            env_mat(:,:,idx) = env;
            idx = idx+1;
        end
    end
    mean_env{n} = mean(env_mat,3);
    clear env_mat idx
end

th_les=[];
r_les=[];

th_bg=[];
r_bg=[];
for r_i = 1:length(r)
    x = r(r_i).*sin(th_scan);
    z = r(r_i).*cos(th_scan);
    
    tmp = find((x.^2+(z-0.0425).^2)< 0.0049^2);
    th_les = [th_les tmp];
    r_les = [r_les r_i.*ones(1,length(tmp))];
    
    tmp1 = find((x.^2+(z-0.0425).^2)< 0.0051^2);
    th_bg = [th_bg tmp1];
    r_bg = [r_bg r_i.*ones(1,length(tmp1))];
    
    clear tmp
end

% mean_env_db = 20*log10(mean_env./max(mean_env(:)));

% figure
% imagesc(mean_env); colormap gray
for n = 1:length(df)
    
    im_les = NaN(size(mean_env{n}));
    im_bg = mean_env{n};

    im_bg([1:1000 3300:end],:) = NaN;
    im_bg(:,[1:20 120:end]) = NaN;
%     im_bg([1:1200 3000:end],:) = NaN;
%     im_bg(:,[1:20 120:end]) = NaN;

    for i = 1:length(r_les)
        im_les(r_les(i),th_les(i)) = mean_env{n}(r_les(i),th_les(i));
    end
    
    for i = 1:length(r_bg)
        im_bg(r_bg(i),th_bg(i)) = NaN;
    end
    
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
%     subplot(121)
%     imagesc(rad2deg(th_scan),100*r,im_les); colormap gray
%     xlabel('\theta (deg)'),ylabel('Depth (cm)')
%     subplot(122)
%     imagesc(rad2deg(th_scan),100*r,im_bg); colormap gray
%     xlabel('\theta (deg)'),ylabel('Depth (cm)')
    
    clear im_les im_bg
end

env_maxCNR = 20*log10(mean_env{3}/max(mean_env{3}(:)));
env_conCNR = 20*log10(env_cont/max(env_cont(:)));

figure
imagesc(rad2deg(th_scan),100*r,env_maxCNR,[-40 0]); colormap gray; axis square
xlabel('Angle (degrees)'),ylabel('Depth (cm)')
% print -dpng -r300 ./hw6_compound.png

figure 
imagesc(rad2deg(th_scan),100*r,env_conCNR,[-40 0]); colormap gray; axis square
xlabel('Angle (degrees)'),ylabel('Depth (cm)')
% print -dpng -r300 ./hw6_control.png

CNR_ratio = CNR(2:end)./CNR(1);

figure
plot(df(2:end),CNR_ratio)
xlabel('\Deltaf (MHz)'),ylabel('Contrast Ratio (CNR_a_v_g/CNR_o_r_i_g)')
% print -dpng -r300 ./hw6_CNR.png

% print -dpng -r300 ./hw6_CNR.png
% figure
% imagesc(im_les); colormap gray
% figure
% imagesc(im_bg); colormap gray

