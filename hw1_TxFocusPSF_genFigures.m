clear all; close all; clc;
label = 'hhp';
load(['sim_data_' label '.mat'])
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

%% PSF B-mode side by side
close all
writerObj1 = VideoWriter(['comparePSF_' label '.avi'])
writerObj1.FrameRate = 10;
open(writerObj1);

figure
for i = 1:length(z_array)
    subplot(121)
    imagesc(x*100,z{1,i}*100,env_db{1,i},[-40 0]);axis image; colormap gray;
%     t1 = title(sprintf('Focused Tx (z = %.2f cm)',100*z_array(i)),'fontsize',16,'Color','w')
    uicontrol('Style', 'text',...
       'String', 'Focused Tx',... %replace something with the text you want
       'Units','normalized',...
       'ForegroundColor','w',...  
       'BackgroundColor','k',...
       'Fontsize',18,...
       'Position', [0.2 0.8 0.2 0.07]); 
   
%     set(t1,'position',[250 300 0])
    axis off;
    
    uicontrol('Style', 'text',...
       'String',  sprintf('z = %.2f cm',100*z_array(i)),... %replace something with the text you want
       'Units','normalized',...
       'ForegroundColor','w',...  
       'BackgroundColor','k',...
       'Fontsize',24,...
       'Position', [0.4 0.25 0.2 0.07]); 
    
    subplot(122)
    imagesc(x*100,z{2,i}*100,env_db{2,i},[-40 0]);axis image; colormap gray;
%     t2 = title(,'fontsize',16,'Color','w')
    uicontrol('Style', 'text',...
       'String', 'Unfocused Tx',... %replace something with the text you want
       'Units','normalized',...
       'ForegroundColor','w',...  
       'BackgroundColor','k',...
       'Fontsize',18,...
       'Position', [0.65 0.8 0.2 0.07]); 
    set(gcf,'color','k');
%     set(t2,'position',[750 300 0])
    axis off
    set(gcf,'position',[100 100 1000 400])
    frame = getframe(gcf);
    writeVideo(writerObj1,frame);
end
close(writerObj1)

%% Resolution Contours
writerObj3 = VideoWriter(['resContours_' label '.avi']);
writerObj3.FrameRate = 5;
writerObj3
open(writerObj3);


for i = 1:length(z_array)
    figure(3)
    clf
    subplot(211)
    hold on
    contour(x*100,z{1,i}*100,env_db{1,i},[-6 -6],'g','linewidth',2);axis image;
    contour(x*100,z{1,i}*100,env_db{1,i},[-12 -12],'m','linewidth',2);axis image;
    ylim([z_array(i)*100 z_array(i)*100+0.15])
    xlim([-0.35 0.35])
    set(gca,'YTickLabel',[]);
    grid on
    xlabel('x (cm)','Fontsize',16)
    ylabel(sprintf('z = %.2f cm',100*z_array(i)),'Fontsize',16)
    title('Focused Tx','Fontsize',18)
    hold off
    subplot(212)
    hold on
    contour(x*100,z{2,i}*100,env_db{2,i},[-6 -6],'g','linewidth',2);axis image;
    contour(x*100,z{2,i}*100,env_db{2,i},[-12 -12],'m','linewidth',2);axis image;
    ylim([z_array(i)*100 z_array(i)*100+0.15])
    xlim([-0.35 0.35])
    set(gca,'YTickLabel',[]);
    grid on
    xlabel('x (cm)','Fontsize',16)
    ylabel(sprintf('z = %.2f cm',100*z_array(i)),'Fontsize',16)
    title('Unfocused Tx','Fontsize',18)
    hold off
    set(gcf,'color','w');
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
xlabel('Depth (cm)','Fontsize',16)
ylabel('Lateral Resolution (cm)','Fontsize',16)
xlim([1.5 8])
legend([p1 p2 p3 p4],{'Focused Tx (6 dB)','Focused Tx (12 dB)','Unfocused Tx (6 dB)','Unfocused Tx (12 dB)'},'Fontsize',14)




% 
% figure
% c1 = contour(x,z,env_db,[-3 -6]);
% grid on 
% % clabel(c1)