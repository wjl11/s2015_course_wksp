clear all; close all; clc;
load ../mat_files/p_data_tissue.mat

% for i=1:5:size(p2,1)
%   imagesc(squeeze(p2(i,:,:))', [-1 1]*p0), title(sprintf('t = %.1f us ',1e6*t(i))), drawnow
%   pause(0.1)
% end

% iline  = round(size(p2,2)/2);
% for i = 1:size(p2,1)
%    	plot(squeeze(p2(i,iline,:))), ylim([-1 1]*p0); title(sprintf('t = %.1f us ',1e6*t(i))), drawnow
% end

c = 1540;
for i = 500
    tmp = squeeze(p2(i,round(size(p2,2)/2),:));
    Fs = 1/dZ
    f = (c)*Fs/2*linspace(0,1,length(tmp)/2+1);
    mag = abs(fft(tmp))/length(tmp);
   	plot(f,mag(1:length(tmp)/2+1)), title(sprintf('t = %.1f us ',1e6*t(i))), drawnow
    [a b] = butter([])

end