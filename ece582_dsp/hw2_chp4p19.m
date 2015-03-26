clear all; close all; clc;

%% 4.19 part a
d0 = 1;
a1 = -1.6454;
a2 = 0.9025;

% Generate realization governed by the autoregressive model:
% AR(2): x(n) = -a1*x(n-1)-a2*x(n-2)+w(n)
x = zeros(1,101);
ni  = 0:100;
for n = ni;
    x(n+1) = 0;
    if n > 1
        x(n+1) = -a1*x(n+1-1);
    end
    if n > 2
        x(n+1) = x(n+1)-a2*x(n+1-2);
    end
    x(n+1) = x(n+1) + randn;
end
figure   
stem(ni,x); xlabel('n'), ylabel('x[n]')

li = 0:100; % calculate the estimated normalized autocorrelation (eqn 1.2.1)
for l = li
    rho_hat(l+1) = 0;
    for n = l:max(ni)
        rho_hat(l+1) = rho_hat(l+1) + x(n+1)*conj(x(n-l+1));
    end
    rho_hat(l+1) = rho_hat(l+1)/sum(x.*conj(x));
end

figure
subplot(211)
stem(li,rho_hat); 
xlabel('$$l$$','Interpreter','Latex')
ylabel('$$\hat{\rho_x}(l)$$','Interpreter','Latex')
xlim([0 10])
ylim([-1 1])
        
%% 4.19 part b
p = roots([1 a1 a2]); % calculate poles of H(z) given a1 and a2
r = unique(abs(p));
theta = unique(abs(angle(p))); % extract angle and magnitude of complex poles

% calculate theoretical value for rho (eqn. 4.2.83)
rho = r.^li.*(sin((li+1).*theta)-r^2*sin((li-1).*theta))./((1+r^2)*sin(theta));
subplot(212)
stem(li,rho,'r'); 
xlabel('$$l$$','Interpreter','Latex')
ylabel('$$\rho_x(l)$$','Interpreter','Latex')
xlim([0 10])
ylim([-1 1])

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
P(1,:) = rho(1:end-1);
for i = 1:99
    P(i+1,:) = [fliplr(rho(2:i+1)) rho(1:100-i)];
end
a_vec = -inv(P)*rho_vec;

figure
hold on
stem(ahat_vec); 
stem(a_vec,'r:')
xlabel('k'); ylabel('Estimated a_k')
hold off
legend('Yule-Walker estimated parameters','True parameters')
xlim([0 15])




