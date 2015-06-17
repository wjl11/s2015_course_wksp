function [b,a] = ldr2dir(K,V)
% Lattice/Ladder form to IIR Direct form Conversion
% -------------------------------------------------
% [b,a] = ldr2dir(K,V)
%  b = numerator polynomial coefficients
%  a = denominator polymonial coefficients
%  K = Lattice coefficients (reflection coefficients)
%  V = Ladder coefficients
%
%-----------------------------------------------------------
% Copyright 2000, by Dimitris G. Manolakis, Vinay K. Ingle,
% and Stephen M. Kogon.  For use with the book
% "Statistical and Adaptive Signal Processing"
% McGraw-Hill Higher Education.
%-----------------------------------------------------------

N = length(K); M = length(V);
V = [V, zeros(1,N-M+1)];
J = 1; a = 1; A = zeros(N,N);
for m=1:1:N
     a = [a,0]+conv([0,K(m)],J);
     A(m,1:m) = -a(2:m+1);
     J = fliplr(a);
end
b(N+1) = V(N+1);
for m = N:-1:1
     A(m,1:m) = A(m,1:m)*V(m+1);
     b(m) = V(m) - sum(diag(A(m:N,1:N-m+1)));
end
