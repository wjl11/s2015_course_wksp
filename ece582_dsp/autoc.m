function r=autoc(x,L)
%
% Computation of the  autocorellation 
% sequence r(l), l=0,1,...,L-1, L <= Lx.
% The mean value is removed
%
% Programmed by: Dimitris Manolakis, 10/4/93
%
%-----------------------------------------------------------
% Copyright 2000, by Dimitris G. Manolakis, Vinay K. Ingle,
% and Stephen M. Kogon.  For use with the book
% "Statistical and Adaptive Signal Processing"
% McGraw-Hill Higher Education.
%-----------------------------------------------------------


Lx=length(x);
x1=zeros(Lx+L-1,1);
x2=x1;
x=x-mean(x);
x1(1:Lx,1)=x;
for k=1:L
	x2=zeros(Lx+L-1,1);
	x2(k:Lx+k-1,1)=x;
	r(k)=x1'*x2;
end
r=r/Lx;
r=r(:);
