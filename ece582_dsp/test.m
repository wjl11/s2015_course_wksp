close all; clear all

rng(5)
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
% 4.2.43. Note the Yule-Walker eqn is modified to accomodate the PZ system
% as is shown in eqn 4.4.9.
khat(1) = -rho_hat(3)/rho_hat(2);
for m = 2:20
    rhohat_vec = rho_hat(3:3+m-1)';
    P = length(rhohat_vec);
    Q = 1; 
    r1 = [fliplr(rho_hat(1:Q+1)) rho_hat(Q+1:Q+1+m-3)];
    c1 = rho_hat(Q+1:Q+1+m-1);
    Phat = toeplitz(c1,r1);
    tmp = -inv(Phat)*rhohat_vec;
    khat(m) = tmp(end);
end