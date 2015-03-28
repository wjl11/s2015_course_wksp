clear all; close all; clc;
rng(0)
%% 4.3 part b and c
x = zeros(1,1000);
X = linspace(-5,5,100);
varx = 2;
ni = 0:999;
for n = ni
    if n >=2
    x(n+1) = x(n+1)+(1/4)*x(n+1-2);
    end
    x(n+1) = x(n+1)+sqrt(15/8)*randn;
end
Y = normpdf(X,0,sqrt(varx));
figure
[counts, centers] = hist(x,50);
hist(x,50)
xlabel('x'), ylabel('Frequency'), title('Observed x[n] N = 1000')


figure
hold on
plot(centers,counts./max(counts),'x')
plot(X,Y./max(Y))
hold off
xlabel('x'), ylabel('Normalized PDF f_x(x)')
legend('Estimated','True')


li = 0:20;
rho_hat = calcAC(x,21);
% *** routine for calcAC defined below used to implement eqn 1.2.1
% for l = li
%     rho_hat(l+1) = 0;
%     for n = l:max(ni)
%         rho_hat(l+1) = rho_hat(l+1) + x(n+1)*conj(x(n-l+1));
%     end
%     rho_hat(l+1) = rho_hat(l+1)/sum(x.*conj(x));
% end

% calculate the normalized autocorrelation given by the expression
r_0 = (1/2).^abs(0)+(-1/2).^abs(0); 
rho = ((1/2).^abs(li)+(-1/2).^abs(li))./r_0; 

figure
hold on
stem(li, rho_hat)
stem(li, rho, 'r:')
hold off
xlabel('$$l$$','Interpreter','Latex')
ylabel('$$\rho_x(l)$$','Interpreter','Latex')
legend('Estimated','True')

clear li rho_hat x rho
%%
pause
close all
%% 4.19 part a
d0 = 1;
a1 = -1.6454;
a2 = 0.9025;

% Generate realization governed by the autoregressive model:
% AR(2): x(n) = -a1*x(n-1)-a2*x(n-2)+w(n)
x = zeros(1,101);
ni  = 0:100;
for n = ni;
    if n >= 1
        x(n+1) = -a1*x(n+1-1);
    end
    if n >= 2
        x(n+1) = x(n+1)-a2*x(n+1-2);
    end
    x(n+1) = x(n+1) + randn;
end
figure   
plot(ni,x); xlabel('n'), ylabel('x[n]')

li = 0:100; % calculate the estimated normalized autocorrelation (eqn 1.2.1)
rho_hat = calcAC(x,101);

figure
hold on
stem(li,rho_hat); 
    
%% 4.19 part b
p = roots([1 a1 a2]); % calculate poles of H(z) given a1 and a2
r = unique(abs(p));
theta = unique(abs(angle(p))); % extract angle and magnitude of complex poles

% calculate theoretical value for rho (eqn. 4.2.83)
rho = r.^li.*(sin((li+1).*theta)-r^2*sin((li-1).*theta))./((1+r^2)*sin(theta));

stem(li,rho,'r:'); 
xlabel('$$l$$','Interpreter','Latex')
ylabel('$$\rho_x(l)$$','Interpreter','Latex')
xlim([0 10])
ylim([-1 1])
legend('Estimated','True')
hold off
%%
pause
close all
%% 4.19 part c
% given the relationships derived from the Yule-Walker Eqns (eqn 4.2.33
% etc)

rhohat_vec = rho_hat(2:end)';
Phat(1,:) = rho_hat(1:end-1);
for i = 1:99
    Phat(i+1,:) = [fliplr(rho_hat(2:i+1)) rho_hat(1:100-i)];
end
ahat_vec = -inv(Phat)*rhohat_vec;

rho_vec = rho(2:end)';
P(1,:) = rho(1:end-1); % form correlation matrix
for i = 1:99
    P(i+1,:) = [fliplr(rho(2:i+1)) rho(1) conj(rho(2:100-i))];
end
a_vec = -inv(P)*rho_vec;

figure
hold on
stem(ahat_vec); 
stem(a_vec,'r:')
xlabel('k'); ylabel('a_k')
hold off
legend('Yule-Walker estimated parameters','True parameters')
xlim([1 15])

%%
pause
close all
%% 4.19 part d

clear tmp

w = -pi:pi/100:pi;
% calculate the PSD given the estimated a_k parameters by definition of Rx
for wi = 1:length(w)
    % Rx = var * |H(z)|^2 where H(z) = 1/(1+a1*z^-1+a2*z^-2 ...)
    tmp = exp(-1i.*w(wi).*(1:length(ahat_vec)));
    R_hat(wi) = 1./((1+sum(ahat_vec'.*tmp))*conj(1+sum(ahat_vec'.*tmp)));
    clear tmp
end

% using eqn 4.2.88, we solve for the true PSD using the poles of the AP system
% with a1 and a2 given by the model
R = d0^2./((1-2.*r.*cos(w-theta)+r^2).*(1-2.*r.*cos(w+theta)+r^2));

figure;
hold on
plot(w,R_hat)
plot(w,R,'r:')
hold off
legend('Estimate PSD','True PSD')
xlim([-pi pi])
xlabel('\omega')
ylabel('PSD')

%%
pause
close all
%% 4.32 part a

clear li ni n l rho_hat x rhohat_vec Phat ahat_vec
ni = 0:99;
x = zeros(1,length(ni));
for n = ni
    if n >= 1
        x(n+1) = x(n+1)+1.9*x(n+1-1);
    end
    if n >= 2
        x(n+1) = x(n+1)-0.9*x(n+1-2);
    end
    x(n+1) = x(n+1) + randn;
end

figure
plot(ni,x); xlabel('n'),ylabel('x[n]')

li = 0:99; % calculate the estimated normalized autocorrelation (eqn 1.2.1)
rho_hat = calcAC(x,length(li));

figure
stem(li, rho_hat); 
xlabel('$$l$$','Interpreter','Latex')
ylabel('$$\hat{\rho_x}(l)$$','Interpreter','Latex')
xlim([0 20])

% calculate PACS based on eqn 4.2.43 
rhohat_vec = rho_hat(2:end)';
Phat(1,:) = rho_hat(1:end-1);
for i = 1:98
    Phat(i+1,:) = [fliplr(rho_hat(2:i+1)) rho_hat(1) conj(rho_hat(2:99-i))];
end
ahat_vec = -inv(Phat)*rhohat_vec; 

figure
hold on
stem(ahat_vec); 
xlabel('m'); ylabel('k_m')
xlim([1 20])

%%
pause
close all
%% 4.32 part b
clear li ni n l rho_hat x rhohat_vec Phat ahat_vec

ni = 0:99;
x = zeros(1,length(ni));

for n = ni
    if n >= 1
        x(n+1) = x(n+1)+1.9*x(n+1-1)-0.5*W(n);
    end
    if n >= 2
        x(n+1) = x(n+1)-0.9*x(n+1-2);
    end
    W(n+1) = randn;
    x(n+1) = x(n+1) + W(n+1);
end

figure
plot(ni,x); xlabel('n'),ylabel('x[n]')

li = 0:99; % calculate the estimated normalized autocorrelation (eqn 1.2.1)
rho_hat = calcAC(x,length(li));

figure
stem(li, rho_hat); 
xlabel('$$l$$','Interpreter','Latex')
ylabel('$$\hat{\rho_x}(l)$$','Interpreter','Latex')
xlim([0 20])

% calculate PACS based on the modified Yule-Walker (eqn 4.4.9) for PZ  
rhohat_vec = rho_hat(3:end)';

% generate the non-hermetian topelitz matrix using the estimated rho values
for i = 1:length(rhohat_vec)-1
    Phat(i,:) = [fliplr(rho_hat(1:i+1)) zeros(1,length(rhohat_vec)-i-1)];
end
Phat(length(rhohat_vec),:) = [fliplr(rho_hat(3:i+2)) rho_hat(2)];
ahat_vec = -inv(Phat)*rhohat_vec;

figure
hold on
stem(ahat_vec); 
xlabel('m'); ylabel('k_m')
xlim([1 20])

%%
pause
%% 5.2 part a & b
close all

N = [11 31 51];
figure
for nn = 1:length(N)
    subplot(length(N),1,nn)
    stem(-floor(N(nn)/2):floor(N(nn)/2),bartlett(N(nn)))
    xlabel('n')
    ylabel('w_B[n]')
    title(sprintf('bartlett (N = %d)',N(nn)))
end

figure
for nn = 1:length(N)
    subplot(length(N),1,nn)
    stem(-floor(N(nn)/2):floor(N(nn)/2),triang(N(nn)))
    xlabel('n')
    ylabel('w_T[n]')
    title(sprintf('triang (N = %d)',N(nn)))
end

Ndtft = 1024;
figure
for nn = 1:length(N)
    subplot(length(N),1,nn)
    WB = fftshift(fft(bartlett(N(nn)),Ndtft));
    WT = fftshift(fft(triang(N(nn)),Ndtft));
    hold on
    plot(linspace(-pi,pi,Ndtft),log10(abs(WB)./max(abs(WB))))
    plot(linspace(-pi,pi,Ndtft),log10(abs(WT)./max(abs(WT))),'r--')
    hold off
    if nn == 1, legend('bartlett','triang'); end;
    xlim([-pi,pi])
    xlabel('\omega (rads)')
    ylabel('Magnitude (dB)')
    title(sprintf('DTFT N = %d',N(nn)))
end

figure
N = [51 81 101 121];
for nn = 1:length(N)
    subplot(length(N),1,nn)
    WR = fftshift(fft(ones(1,51),Ndtft));
    WB = fftshift(fft(bartlett(N(nn)),Ndtft));
    hold on
    plot(linspace(-pi,pi,Ndtft),log10(abs(WR)./max(abs(WR))))
    plot(linspace(-pi,pi,Ndtft),log10(abs(WB)./max(abs(WB))),'r--')
    hold off
    if nn == 1, legend('bartlett','rect'); end;
    xlim([-pi/6,pi/6])
    xlabel('\omega (rads)')
    ylabel('Magnitude (dB)')
    title(sprintf('DTFT N = %d',N(nn)))
end
