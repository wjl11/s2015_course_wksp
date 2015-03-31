function [K] = df2latcf(a)
% FIR Direct form to All-Zero Lattice form Conversion
% ---------------------------------------------------
% [K] = df2latcf(a)
%  K = Lattice filter coefficients (reflection coefficients)
%  a = FIR direct form coefficients (impulse response)
%
%-----------------------------------------------------------
% Copyright 2000, by Dimitris G. Manolakis, Vinay K. Ingle,
% and Stephen M. Kogon.  For use with the book
% "Statistical and Adaptive Signal Processing"
% McGraw-Hill Higher Education.
%-----------------------------------------------------------


M = length(a);
K = zeros(1,M);
a1 = a(1);
if a1 == 0
	error('a(1) is equal to zero')
end
K(1) = a1; A = a/a1;
for m=M:-1:2
	K(m) = A(m);
	J = fliplr(A);
	A = (A-K(m)*J)/(1-K(m)*K(m));
	A = A(1:m-1);
end
