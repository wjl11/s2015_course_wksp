clear all; close all; clc;
load('sim_data.mat')
%% PSF B-mode
% writerObj1 = VideoWriter('focusTxPSF.avi')
writerObj1.FrameRate = 10;
open(writerObj1);

% writerObj2 = VideoWriter('pwTxPSF.avi')
writerObj2.FrameRate = 10;
open(writerObj2);

for i = 1:length(z_array)
    figure(1)
    imagesc(x*100,z{1,i}*100,env_db{1,i},[-40 0]);axis image; colormap gray;
    title(sprintf('Focused Tx (z = %.2f cm)',100*z_array(i)),'Color','w')
    axis off;
    set(gcf,'color','k');
    
    frame = getframe(gcf);
    writeVideo(writerObj1,frame);
end
close(writerObj1)

for i = 1:length(z_array)
    figure(2)
    imagesc(x*100,z{2,i}*100,env_db{2,i},[-40 0]);axis image; colormap gray;
    title(sprintf('Unfocused Tx (z = %.2f cm)',100*z_array(i)),'Color','w')
    axis off;
    set(gcf,'color','k');
    
    frame = getframe(gcf);
    writeVideo(writerObj2,frame);

end
close(writerObj2)


%% Resolution Contours
writerObj3 = VideoWriter('resContours.avi');
writerObj3.FrameRate = 3;
writerObj3
open(writerObj3);

figure(3)
for i = 1:length(z_array)
    clf
    subplot(211)
    hold on
    contour(x*100,z{1,i}*100,env_db{1,i},[-6 -6],'g');axis image;
    contour(x*100,z{1,i}*100,env_db{1,i},[-12 -12],'m');axis image;
    ylim([z_array(i)*100 z_array(i)*100+0.15])
    xlim([-0.35 0.35])
    set(gca,'YTickLabel',[]);
    grid on
    xlabel('x (cm)')
    ylabel(sprintf('z = %.2f cm',100*z_array(i)))
    title('Focused Tx')
    hold off
    subplot(212)
    hold on
    contour(x*100,z{2,i}*100,env_db{2,i},[-6 -6],'g');axis image;
    contour(x*100,z{2,i}*100,env_db{2,i},[-12 -12],'m');axis image;
    ylim([z_array(i)*100 z_array(i)*100+0.15])
    xlim([-0.35 0.35])
    set(gca,'YTickLabel',[]);
    grid on
    xlabel('x (cm)')
    ylabel(sprintf('z = %.2f cm',100*z_array(i)))
    title('Unfocused Tx')
    hold off
    frame = getframe(gcf);
    writeVideo(writerObj3,frame);
%     disp('testing')
end
close(writerObj3)

%% Resolution vs Depth
figure
hold on
p1=plot(z_array*100,fw6db(:,1),'b','linewidth',2);
p2=plot(z_array*100,fw12db(:,1),'b:','linewidth',2);
p3=plot(z_array*100,fw6db(:,2),'r','linewidth',2)
p4=plot(z_array*100,fw12db(:,2),'r:','linewidth',2)
grid on
hold off
xlabel('Depth (cm)')
ylabel('Lateral Resolution (cm)')
xlim([2 8])
legend([p1 p2 p3 p4],{'Focused Tx (6 dB)','Focused Tx (12 dB)','Unfocused Tx (6 dB)','Unfocused Tx (12 dB)'})




% 
% figure
% c1 = contour(x,z,env_db,[-3 -6]);
% grid on 
% % clabel(c1)