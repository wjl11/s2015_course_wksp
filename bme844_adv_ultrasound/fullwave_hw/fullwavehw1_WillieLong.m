%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GIANMARCO PINTON
% WRITTEN: NOV 13, 2013
% LAST MODIFIED: NOV 13, 2013
% Launch Fullwave code, easy matlab wrapper
%%% Basic variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear all; clc;
c0 = 1540;         % speed of sound (m/s)
omega0 = 2*pi*1e6; % center radian frequency of transmitted wave
wY = 2e-2;         % width of simulation field (m)
wZ = 3e-2;         % depth of simulation field (m)
duration = 40e-6;  % duration of simulation (s)
p0 = 1e5; % pressure in Pa
%%% Advanced variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ppw = 15;           % number of points per spatial wavelength
cfl = 0.4;         % Courant-Friedrichs-Levi condition
%%% Grid size calculations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lambda = c0/omega0*2*pi;
nY = round(wY/lambda*ppw);  % number of lateral elements
nZ = round(wZ/lambda*ppw);  % number of depth elements
nT = round(duration*c0/lambda*ppw/cfl);
dY = c0/omega0*2*pi/ppw;
dZ = c0/omega0*2*pi/ppw;
dT = dY/c0*cfl;
%%% Generate field maps %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cmap = ones(nY,nZ)*1540;   % speed of sound map (m/s)
rhomap = ones(nY,nZ)*1000; % density map (kg/m^3)
Amap = ones(nY,nZ)*0.0;    % attenuation map (dB/MHz/cm)
boveramap = -2*ones(nY,nZ);    % nonlinearity map 
cmap(round(nY/3)-1:round(nY/3)+1,round(nZ/1.3)-1:round(nZ/1.3)+1)=0.5*c0; % scatterer (EDITED for prob. 1)
imagesc(cmap'), axis equal, axis tight
%%% Generate input coordinates %%%%%%%%%%%%%%%%%%%%%%%%%%%
inmap = zeros(nY,nZ); 
inmap(:,1) = ones(nY,1); inmap(:,2) = ones(nY,1); inmap(:,3) = ones(nY,1);
imagesc(inmap'), axis equal, axis tight
incoords = mapToCoords(inmap); % note zero indexing for compiled code
plot(incoords(:,1),incoords(:,2),'.')
%%% Generate initial conditions based on input coordinates %%%%%%
ncycles = 2; % number of cycles in pulse
dur = 2; % exponential drop-off of envelope

fcen=[round(nY/3) round(nZ/1.3)]; % center of focus (EDITED for prob. 1)
t = (0:nT-1)/nT*duration-ncycles/omega0*2*pi;
icmat=zeros(size(incoords,1),nT);
icvec = exp(-(1.05*t*omega0/(ncycles*pi)).^(2*dur)).*sin(t*omega0)*p0;
plot(icvec), hold all
icmat(1:size(incoords,1)/3,:) = focusCoords(fcen(1),fcen(2),incoords(1:size(incoords,1)/3,:),icvec,cfl);
t=t-dT/cfl;
icvec = exp(-(1.05*t*omega0/(ncycles*pi)).^(2*dur)).*sin(t*omega0)*p0;
plot(icvec)
icmat(size(incoords,1)/3+1:size(incoords,1)/3*2,:)=focusCoords(fcen(1),fcen(2),incoords(size(incoords,1)/3+1:size(incoords,1)/3*2,:),icvec,cfl);
t=t-dT/cfl;
icvec = exp(-(1.05*t*omega0/(ncycles*pi)).^(2*dur)).*sin(t*omega0)*p0;
plot(icvec), hold off
icmat(size(incoords,1)/3*2+1:size(incoords,1),:)=focusCoords(fcen(1),fcen(2),incoords(size(incoords,1)/3*2+1:size(incoords,1),:),icvec,cfl);
imagesc(icmat)
%%% Generate output coordinates %%%%%%%%%%%%%%%%%%%%%%%%%%
outmap = zeros(nY,nZ);
[modidy modidz] = meshgrid(1:2:nY,1:2:nZ);
outmap(modidy,modidz) = 1;
imagesc(outmap'), axis equal, axis tight
outcoords = mapToCoords(outmap);
%%% Launch %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
launchTotalFullWave2(c0,omega0,wY,wZ,duration,p0,ppw,cfl,cmap',rhomap',Amap',boveramap',incoords,outcoords,icmat);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
!./try6_nomex
toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ncoordsout=size(outcoords,1)  
nRun=sizeOfFile('genout.dat')/4/ncoordsout
genout = readGenoutSlice(['genout.dat'],0:nRun-1,size(outcoords,1));
p = reshape(genout,size(genout,1),size(modidy,2),size(modidy,1));
imagesc(squeeze(p(end,:,:))), colorbar

for i=1:5:size(p,1)
  imagesc(squeeze(p(i,:,:))', [-1 1]*p0), title(sprintf('t = %.1f us ',1e6*t(i))), drawnow
end

i = find(t>=21e-6,1);
y = linspace(-wY/2,wY/2,nY);
z = linspace(0,wZ,nZ);
figure
imagesc(1000.*y,1000.*z,squeeze(p(i,:,:))', [-1 1]*p0), title(sprintf('t = %.1f us ',1e6*t(i))); 
axis equal, axis tight
xlabel('y (mm)'), ylabel('z (mm)')
print -djpeg ./fig1_WillieLong.jpg

%%

