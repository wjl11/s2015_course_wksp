close all; clear all;

r = 3.5;
th = [-20:20:20].*pi/180;
rangex = [-1 1]; % focal range of x [cm]
rangez = [2.5 3.5]; % focal range of z relative to scan rad [cm]
apertureL = [-2 2];
im_depth = [0 5];

% rangex = -1;
% rangez = 3;

if length(rangex) == 1 & length(rangez) ==1
    gridx = rangex;
    gridz = rangez;
else
    gridx = rangex(1):0.5:rangex(2);
    gridz = rangez(1):0.5:rangez(2);
end

focal_pts(:,:,1) = repmat(gridx,[length(gridz) 1]);
focal_pts(:,:,2) = repmat(gridz',[1 length(gridx)]);

x = focal_pts(:,:,1)./100;
z = focal_pts(:,:,2)./100;
r = r/100;

xi_tmp = zeros(size(x));
zi_tmp = zeros(size(z));
xi = zeros(size(x,1), size(x,2), length(th));
zi = zeros(size(z,1), size(z,2), length(th));

focus_rf = zeros(length(gridz),length(gridx),length(th));

for nt = 1:length(th) 
    th_i = atan(x./(r-z));
    % off axis focal points
    xi_tmp(th_i ~= 0) = -(x(th_i ~= 0)./sin(th_i(th_i ~= 0))).*sin(th(nt)-th_i(th_i ~= 0));
    zi_tmp(th_i ~= 0) = r-(x(th_i ~= 0)./sin(th_i(th_i ~= 0))).*cos(th(nt)-th_i(th_i ~= 0));
    % on axis focal points
    xi_tmp(th_i == 0) = -(r-z(th_i == 0)).*sin(th(nt));
    zi_tmp(th_i == 0) = r - (r-z(th_i == 0)).*cos(th(nt));
    % at scan center (z = r and x = 0)
    xi_tmp(isnan(th_i)) = 0;
    zi_tmp(isnan(th_i)) = r;

    xi(:,:,nt) = double(xi_tmp);
    zi(:,:,nt) = double(zi_tmp);
    
    figure;
    plot(xi_tmp(:), zi_tmp(:),'x','Linewidth',2); axis image
    set(gca,'YDir','reverse')
    xlim(apertureL.*1e-2)
    ylim(im_depth.*1e-2)
    title(['Focal points for frame at ' num2str(th(nt)*180/pi) ' degrees'])

end

figure
plot(xi(:),'x')

figure
plot(zi(:),'x')