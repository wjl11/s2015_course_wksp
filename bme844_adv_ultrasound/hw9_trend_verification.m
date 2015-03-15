%%
close all; clear all 
x = -2*pi:pi/1000:2*pi;

y = sinc(x).^2;

det_y = abs(hilbert(y));

figure
plot(x,y,x,det_y)

figure
hold on
plot(fftshift(abs(fft(y))))
plot(fftshift(abs(fft(det_y))))
hold off

%%
close all; clear all;
load('./mat_files/hw2_foc_psf_v2.mat');
k_rf = fftshift(abs(fft2(rf{1})));
k_env = fftshift(abs(fft2(env{1})));
figure
subplot(211)
imagesc(k_rf)
title('RF k-space')
subplot(212)
imagesc(k_env)
title('Detected k-space')

figure
subplot(211)
hold on
plot(k_rf(310,:))
plot(k_env(340,:))
hold off
subplot(212)
hold on
p1 = plot(k_rf(:,85))
p2 = plot(k_env(:,85))
hold off
legend([p1 p2],{'RF','Detected'})

%% speckle k-space
close all; clear all; 
load('./mat_files/hw5_data_full.mat');
k_rf = fftshift(abs(fft2(rf{10,10}(300:600,:))));
k_env = fftshift(abs(fft2(env{10,10}(300:600,:))));

figure
subplot(211)
imagesc(k_rf)
title('RF k-space')
subplot(212)
imagesc(k_env)
title('Detected k-space')

figure
subplot(211)
hold on
plot(k_rf(134,:))
plot(k_env(151,:))
hold off
subplot(212)
hold on
p1 = plot(k_rf(:,53))
p2 = plot(k_env(:,53))
hold off
legend([p1 p2],{'RF','Detected'})