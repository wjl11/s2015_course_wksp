clear all; close all; clc;
load ../mat_files/p_data_tissue.mat

test = 200;
figure
for i=test
  imagesc(squeeze(p2(i,:,:))', [-1 1]*p0), title(sprintf('t = %.1f us ',1e6*t(i))), drawnow
%   pause(0.1)
end

var = 0.05e6;
c = 1540;
fsZ = 1/dZ;
fsY = 1/dY;
fz = c*fsZ/2*linspace(-1,1,size(p2,3));
fy = c*fsY/2*linspace(-1,1,size(p2,2));
filt_f0 = gaussmf(fz,[var 1e6]);
filt_f1 = gaussmf(fz,[var 2e6]);

figure
% for i=1:size(p2,1)
for i=1:10:size(p2,1)
    tmp = fft2(squeeze(p2(i,:,:))');
    mag = fftshift(abs(tmp));
%     imagesc(fy,fz,mag), drawnow, pause(0.1);
    plot(fz,mag(:,find(fy == 0))),ylim([0 1e8]), drawnow, pause(0.2); 
%     tmp = squeeze(p2(i,round(size(p2,2)/2),:));
%     mag = abs(fft(tmp))/length(tmp);
%     mag = mag(1:length(tmp)/2+1)';
   	
% 
%     E_f0(i) = sum(mag.*filt_f0);
%     E_f1(i) = sum(mag.*filt_f1);
%     
%     plot(f,mag.*filt_f1), title(sprintf('t = %.1f us ',1e6*t(i))), drawnow
end
% 
% hold on
% plot(t,E_f0./max(E_f0))
% plot(t,E_f1./max(E_f1))
% hold off