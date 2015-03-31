%% 5.2 part a & b
clear all; clc; close all

% generate and plot bartlett windows of varying length N
N = [11 31 51];
figure
for nn = 1:length(N)
    subplot(length(N),1,nn)
    stem(-floor(N(nn)/2):floor(N(nn)/2),bartlett(N(nn)))
    xlabel('n')
    ylabel('w_B[n]')
    title(sprintf('bartlett (N = %d)',N(nn)))
end

% generate and plot triang windows of varying length N
figure
for nn = 1:length(N)
    subplot(length(N),1,nn)
    stem(-floor(N(nn)/2):floor(N(nn)/2),triang(N(nn)))
    xlabel('n')
    ylabel('w_T[n]')
    title(sprintf('triang (N = %d)',N(nn)))
end

% generate plots of dtft of bartlett and triang windows of varying N
Ndtft = 4*1024;
figure
for nn = 1:length(N)
    % calculate the dtft of the bartlett windows in dB with zero pad
    tmp = fft(bartlett(N(nn)),Ndtft);
    WB(nn,:) = 10.*log10(abs(tmp)./max(abs(tmp)));
    clear tmp
    
    % calculate the dtft of the triang windows in dB with zero pad
    tmp = fft(triang(N(nn)),Ndtft);
    WT(nn,:) = 10.*log10(abs(tmp)./max(abs(tmp)));
    clear tmp
    
    % plot the dtft of windows in dB
    subplot(length(N),1,nn)
    hold on
    plot(linspace(-pi,pi,Ndtft),fftshift(WB(nn,:)));
    plot(linspace(-pi,pi,Ndtft),fftshift(WT(nn,:)),'r--');
    hold off
    if nn == 1, legend('bartlett','triang'); end;
    xlim([-pi,pi])
    xlabel('\omega (rads)')
    ylabel('Magnitude (dB)')
    title(sprintf('DTFT N = %d',N(nn)))
end

% generate plot of main lobe width as a function of window length
N = 11:5:81
for nn = 1:length(N)
    % calculate dtfts of windows for varying N
    tmp = fft(bartlett(N(nn)),Ndtft);
    WB(nn,:) = 10.*log10(abs(tmp)./max(abs(tmp)));
    clear tmp
    
    tmp = fft(triang(N(nn)),Ndtft);
    WT(nn,:) = 10.*log10(abs(tmp)./max(abs(tmp)));
    clear tmp
    w = linspace(0,2*pi,Ndtft);
    
    % find the frequency value corresponding to full width half max of
    % bartlett and triang windows of varying N (where db = -3)
    fwhmB(nn) = 2*w(find(WB(nn,:)<=-3,1));
    fwhmT(nn) = 2*w(find(WT(nn,:)<=-3,1));
    
    % find the frequency value corresponding to null to null main lobe
    % width for bartlett and triang windows of varying N
    [~, itmp] = findpeaks(-WB(nn,:));
    bwB(nn) = 2*w(itmp(1));
    [~, itmp] = findpeaks(-WT(nn,:));
    bwT(nn) =  2*w(itmp(1));
end
figure; hold on
plot(N,bwB,'x')
plot(N,bwT,'rx')
hold off
xlabel('Window Length N')
ylabel('Null-to-null BW (rads)')
title('Mainlobe Width for Varying N')
legend('bartlett','triang')

%% 5.2 part c

clear WB bwB
% calculate and plot the mainlobe width for series of bartlett windows of
% varying N and compare to mainlobe width of 51-point rect
N = 51:5:121
for nn = 1:length(N)
    % calculate dtfts of windows for varying N
    tmp = fft(bartlett(N(nn)),Ndtft);
    WB(nn,:) = 10.*log10(abs(tmp)./max(abs(tmp)));
    clear tmp
    
    tmp = fft(ones(1,51),Ndtft);
    WR(nn,:) = 10.*log10(abs(tmp)./max(abs(tmp)));
    clear tmp
    w = linspace(0,2*pi,Ndtft);

    % find the frequency value corresponding to null to null main lobe
    % width for bartlett and rect windows of varying N
    [~, itmp] = findpeaks(-WB(nn,:));
    bwB(nn) = 2*w(itmp(1));
    [~, itmp] = findpeaks(-WR(nn,:));
    bwR(nn) =  2*w(itmp(1));
end
figure; hold on
plot(N,bwB,'x')
plot(N,bwR,'r--')
hold off
xlabel('Window Length N')
ylabel('Null-to-null BW (rads)')
title('Mainlobe Width for Varying N')
legend('bartlett','rect')

clear WB bwB
% compare dtft of bartlett of varying lengths to rect window with N = 51
figure
N = [51 81 101 121];
for nn = 1:length(N)
    
    % calculate dtft of rect and bartlett windows
    WR = fft(ones(1,51),Ndtft);
    WB = fft(bartlett(N(nn)),Ndtft);
    
    subplot(length(N),1,nn)
    hold on
    plot(linspace(-pi,pi,Ndtft),fftshift(10.*log10(abs(WB)./max(abs(WB)))))
    plot(linspace(-pi,pi,Ndtft),fftshift(10.*log10(abs(WR)./max(abs(WR)))),'r--')
    hold off
    
    if nn == 1, legend('bartlett','rect'); end;
    xlim([-pi/6,pi/6])
    xlabel('\omega (rads)')
    ylabel('Magnitude (dB)')
    title(sprintf('DTFT N = %d',N(nn)))
end


