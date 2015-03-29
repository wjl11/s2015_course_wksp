function rho_hat = calcAC(x,nlags)
% calculate the estimated normalized autocorrelation (eqn 1.2.1)
% nlags specifics number of lags to calculate rho (starting from l = 0)
% x is the signal to autocorr whose 1st index corresponds to n = 0

rho_hat = zeros(1,nlags);
for l = 0:nlags-1
    rho_hat(l+1) = 0;
    for n = l:(length(x)-1)
        rho_hat(l+1) = rho_hat(l+1) + x(n+1)*conj(x(n-l+1));
    end
    rho_hat(l+1) = rho_hat(l+1)/sum(x.*conj(x));
end