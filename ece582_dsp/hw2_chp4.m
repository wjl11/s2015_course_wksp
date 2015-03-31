clear all; close all; clc;
rng(0)
%% 4.3 part b and c

% generate 1000 realizations of x[n] using difference equation
x = zeros(1,1000);
ni = 0:999;
for n = ni
    if n >=2
    x(n+1) = x(n+1)+(1/4)*x(n+1-2);
    end
    x(n+1) = x(n+1)+sqrt(15/8)*randn;
end

% organize x[n] realizations in histogram to estimate pdf
figure
[counts, centers] = hist(x,50);
hist(x,50)
xlabel('x'), ylabel('Frequency'), title('Histogram of x[n] N = 1000')

% generate true pdf values based on calculate mean and variance of x
varx = 2;
mux = 0;
X = linspace(min(centers),max(centers),100);
Y = normpdf(X,mux,sqrt(varx));

% plot normalized estimated and true pdfs
tmp = diff(centers);
dx = tmp(1);
figure
hold on
% estimated histogram is normalized to ensure AUC = 1
plot(centers,counts./(dx*sum(counts)),'x') 
plot(X,Y)
hold off
xlabel('x'), ylabel('Probability f_x(x)')
legend('Estimated','True')
title('Estimated and True PDF')

% estimate the normalized autocorrelation using user defined function
% calcAC, implementing eqn 1.2.1 (see section below)
li = 0:20;
rho_hat = calcAC(x,21);

% calculate the normalized autocorrelation given by the known expression
r_0 = (1/2).^abs(0)+(-1/2).^abs(0); 
rho = ((1/2).^abs(li)+(-1/2).^abs(li))./r_0; 

% plot estimated and true rho values
figure
hold on
stem(li, rho_hat)
stem(li, rho, 'r:')
hold off
xlabel('$$l$$','Interpreter','Latex')
ylabel('$$\rho_x(l)$$','Interpreter','Latex')
legend('Estimated','True')
title('Estimated and True Autocorrelation')

clear li rho_hat x rho
%%
pause
close all
%% 4.19 part a & b

% Generate 100 samples of x[n] given by the autoregressive model:
% AR(2): x(n) = -a1*x(n-1)-a2*x(n-2)+w(n)
d0 = 1;
a1 = -1.6454;
a2 = 0.9025;

x = zeros(1,100);
ni  = 0:99;
for n = ni;
    if n >= 1
        x(n+1) = -a1*x(n+1-1);
    end
    if n >= 2
        x(n+1) = x(n+1)-a2*x(n+1-2);
    end
    w(n+1) = randn;
    x(n+1) = x(n+1) + w(n+1);
end
figure
plot(ni,x); xlabel('n'), ylabel('x[n]')
title('x[n] for N = 100')

% estimate ACS (eqn 1.2.1)
li = 0:99; 
rho_hat = calcAC(x,100);

% calculate poles of H(z) given a1 and a2
p = roots([1 a1 a2]); 

% extract angle and magnitude of complex poles to solve for true rho values
r = unique(abs(p));
theta = unique(abs(angle(p))); 

% calculate theoretical ACS (eqn. 4.2.83)
rho = r.^li.*(sin((li+1).*theta)-r^2*sin((li-1).*theta))...
    ./((1+r^2)*sin(theta));

% plot estimated and theoretical ACS
figure
hold on
stem(li,rho_hat); 
stem(li,rho,'r:'); 
xlabel('$$l$$','Interpreter','Latex')
ylabel('$$\rho_x(l)$$','Interpreter','Latex')
xlim([0 20])
ylim([-1 1])
legend('Estimated','True')
title('Estimated and Theoretical ACS')
hold off
%%
pause
close all
%% 4.19 part c

% define estimated rho vector for Yule-Walker 
rhohat_vec = rho_hat(2:end)';

% define estimated autocorr matrix for Yule-Walker
r1 = conj(rho_hat(1:end-1)); % first row of toeplitz hermitian matrix
Phat = toeplitz(r1); clear r1

% solve for estimated parameters using the Yule-Walker eqn (eqn 4.2.33)
ahat_vec = -inv(Phat)*rhohat_vec;

% define true rho vector for Yule-Walker 
rho_vec = rho(2:end)';

% define true autocorr matrix for Yule-Walker
r1 = conj(rho(1:end-1)); 
P = toeplitz(r1); clear r1

% solve for true parameters using the Yule-Walker eqn (eqn 4.2.33)
a_vec = -inv(P)*rho_vec;

% plot true and estimated parameter values for comparison
figure
hold on
stem(ahat_vec); 
stem(a_vec,'r:')
xlabel('k'); ylabel('a_k')
hold off
legend('Estimated','True')
xlim([0 15])
title('Estimated and Theoretical a_k Parameters')

%%
pause
close all
%% 4.19 part d

clear tmp

w = -pi:pi/100:pi;
% calculate the PSD given the estimated a_k parameters by definition of Rx
% Rx = w_var * |H(z)|^2 where H(z) = 1/(1+a1*z^-1+a2*z^-2 ...) 
%   z = e^jw and w_var = 1 (p. 165)
for wi = 1:length(w)
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
legend('Estimated','True')
xlim([-pi pi])
xlabel('\omega')
ylabel('PSD')
title('Estimated and Theoretical PSD')
%%
pause
close all 
%% 4.19 part e

% calculate the true and estimated PACS/lattice structures from the rho 
% values through application of the Yule-Walker equation given by eqn.
% 4.2.43

for m = 1:15
    % calculate PACS using estimated rho
    % evaluate Yule-Walker using a portion of autocorrelation values where
    % m is the maximum index 
    rhohat_vec = rho_hat(2:m+1)';
    r1 = conj(rho_hat(1:m));
    Phat = toeplitz(r1);
    tmp = -inv(Phat)*rhohat_vec;
    
    % extract the mth value which defines the lattice/PACS parameter
    khat(m) = tmp(end);
    clear r1 c1 tmp
    
    % repeat procedure using true rho values
    rho_vec = rho(2:m+1)';
    r1 = conj(rho(1:m));
    P = toeplitz(r1);
    tmp = -inv(P)*rho_vec;
    k(m) = tmp(end);
    clear r1 c1 tmp
end

figure
hold on
stem(khat); 
stem(k,'r:')
xlabel('m'); ylabel('a^(^m^)_m')
hold off
legend('Estimated','True')
xlim([0 15])
title('Estimated and Theoretical PACS')
%%
pause
%% 4.32 part a
close all; clear all;

% Generate 100 samples of x[n] in part a) based on difference equation 
% defined by the system function H(z)
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
title('x[n] for N = 100')

% estimate the normalized autocorrelation using user defined function
% calcAC, implementing eqn 1.2.1
li = 0:99; 
rho_hat = calcAC(x,length(li));

figure
subplot(211)
stem(li, rho_hat); 
xlabel('$$l$$','Interpreter','Latex')
ylabel('$$\hat{\rho_x}(l)$$','Interpreter','Latex')
xlim([0 20])
title('Estimated ACS')

% calculate the estimated PACS/lattice structures from the rho 
% values through application of the Yule-Walker equation given by eqn.
% 4.2.43
for m = 1:20
    rhohat_vec = rho_hat(2:m+1)';
    r1 = conj(rho_hat(1:m));
    Phat = toeplitz(r1);
    tmp = -inv(Phat)*rhohat_vec;
    khat(m) = tmp(end);
end

subplot(212)
stem(khat); 
xlabel('m'); ylabel('a^(^m^)_m')
xlim([0 20])
title('Estimated PACS')

%%
pause
%% 4.32 part b
close all; clear all

% Generate 100 samples of x[n] in part b) based on difference equation 
% defined by the system function H(z)
ni = 0:99;
x = zeros(1,length(ni));
for n = ni
    if n >= 1
        x(n+1) = x(n+1)+1.9*x(n+1-1)-0.5*W(n+1-1);
    end
    if n >= 2
        x(n+1) = x(n+1)-0.9*x(n+1-2);
    end
    W(n+1) = randn; % note that the system is PZ
    x(n+1) = x(n+1) + W(n+1);
end
figure
plot(ni,x); xlabel('n'),ylabel('x[n]')
title('x[n] for N = 100')

% estimate the normalized autocorrelation using user defined function
% calcAC, implementing eqn 1.2.1
li = 0:99;
rho_hat = calcAC(x,length(li));

figure
subplot(211)
stem(li, rho_hat); 
xlabel('$$l$$','Interpreter','Latex')
ylabel('$$\hat{\rho_x}(l)$$','Interpreter','Latex')
xlim([0 20])
title('Estimated ACS')


% calculate the estimated PACS/lattice structures from the rho 
% values through application of the Yule-Walker equation given by eqn.
% 4.2.43.
for m = 1:20
    rhohat_vec = rho_hat(2:m+1)';
    r1 = conj(rho_hat(1:m));
    Phat = toeplitz(r1);
    tmp = -inv(Phat)*rhohat_vec;
    khat(m) = tmp(end);
end

subplot(212)
stem(khat); 
xlabel('m'); ylabel('a^(^m^)_m')
xlim([0 20])
title('Estimated PACS')