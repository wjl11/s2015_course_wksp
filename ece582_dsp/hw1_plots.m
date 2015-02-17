%% 2.2c
clear all; close all; clc
w0 = pi/4;
N = 32;

A = 1;

w = (-pi:pi/10000:pi);
    
Xr = (A/2).*(cos((w-w0).*(N-1)/2).*(sin((w-w0).*N/2)./sin((w-w0)./2))) ...
    +(A/2).*(cos((w+w0).*(N-1)/2).*(sin((w+w0).*N/2)./sin((w+w0)./2)));
Xi = -(A/2).*(sin((w-w0).*(N-1)/2).*(sin((w-w0).*N/2)./sin((w-w0)./2))) ...
    -(A/2).*(sin((w+w0).*(N-1)/2).*(sin((w+w0).*N/2)./sin((w+w0)./2)));

n = 1:32;
N = length(n);
x_dft = A.*cos(w0*(n-1));
figure
plot(x_dft)
X_dft = fftshift(fft(x_dft));
w_dft = (0:2/N:2*(N-1)/N).*pi-pi;

figure;
hold on
plot(w,Xr)
stem(w_dft,real(X_dft),'r')
title('N = 32, \omega_0 = \pi/4'), ylabel('Real X_R'),xlabel('Frequency \omega')
hold off
% print -dpng -r300 ./hw1_2-2c_real.png

figure
hold on
plot(w,Xi)
stem(w_dft,imag(X_dft),'r')
title('N = 32, \omega_0 = \pi/4'), ylabel('Imaginary X_I'),xlabel('Frequency \omega')
hold off
% print -dpng -r300 ./hw1_2-2c_im.png

%% 2.2d
clear all; close all; clc
w0 = 1.1*pi/4;
N = 32;

A = 1;

w = (-pi:pi/10000:pi);
    
Xr = (A/2).*(cos((w-w0).*(N-1)/2).*(sin((w-w0).*N/2)./sin((w-w0)./2))) ...
    +(A/2).*(cos((w+w0).*(N-1)/2).*(sin((w+w0).*N/2)./sin((w+w0)./2)));
Xi = -(A/2).*(sin((w-w0).*(N-1)/2).*(sin((w-w0).*N/2)./sin((w-w0)./2))) ...
    -(A/2).*(sin((w+w0).*(N-1)/2).*(sin((w+w0).*N/2)./sin((w+w0)./2)));

n = 1:32;
N = length(n);
x_dft = A.*cos(w0*(n-1));
figure
stem(x_dft)
X_dft = fftshift(fft(x_dft));
w_dft = (0:2/N:2*(N-1)/N).*pi-pi;

figure;
hold on
plot(w,Xr)
stem(w_dft,real(X_dft),'r')
title('N = 32, \omega_0 = 1.1\pi/4'), ylabel('Real X_R'),xlabel('Frequency \omega')
hold off

% print -dpng -r300 ./hw1_2-2d_real.png

figure
hold on
plot(w,Xi)
stem(w_dft,imag(X_dft),'r')
title('N = 32, \omega_0 = 1.1\pi/4'), ylabel('Imaginary X_I'),xlabel('Frequency \omega')
hold off

% print -dpng -r300 ./hw1_2-2d_im.png

%% 2.5a
clear all; close all; clc

% analytical solution
n_an = 0:100;
xx_an = (n_an+1).*(0.9).^n_an;

figure
stem(n_an,xx_an)
xlabel('n'), ylabel('x[n]*x[n]'),title('Analytical solution')
% print -dpng -r300 ./hw1_2-5a.png

% convolution via matlab
n = 0:50;
x = (0.9).^n;
xx_conv = conv(x,x);
n_conv = 0:length(xx_conv)-1;
figure
stem(n_conv,xx_conv)
xlabel('n'), ylabel('x[n]*x[n]'),title('Solution via conv()')
% print -dpng -r300 ./hw1_2-5b.png

% filter via matlab
n_filt = 0:100;
x = (0.9).^n_filt;
b = 1;
a = [1 -0.9];
xx_filt = filter(b,a,x);
figure
stem(n_filt,xx_filt)
xlabel('n'), ylabel('x[n]*x[n]'),title('Solution via filter()')
% print -dpng -r300 ./hw1_2-5c.png

figure
hold on
s1 = stem(n_an,abs(xx_an-xx_conv));
s2 = stem(n_an,abs(xx_an-xx_filt),'r');
xlabel('n'), ylabel('Absolute Error x[n]*x[n]')
title('Matlab error (\Deltax[n]*x[n] from analytical solution)')
hold off
legend([s1(1),s2(1)],{'Error from conv()','Error from filter()'})
% print -dpng -r300 ./hw1_2-5d.png

%% 3.7
clear all; close all; clc

y = 0:0.1:10;

fy2 = y.*exp(-y);
fy3 = (y.^2/2).*exp(-y);
fy4 = (y.^3/6).*exp(-y);

sum(fy4)
sum(fy3)

hold on;
p1 = plot(y,fy2)
p2 = plot(y,fy3)
p3 = plot(y,fy4)
hold off

legend([p1 p2 p3],{'f_y_2','f_y_3','f_y_4'})
xlabel('y'), ylabel('PDF f(y)')

% print -dpng -r300 ./hw1_3.7.png
