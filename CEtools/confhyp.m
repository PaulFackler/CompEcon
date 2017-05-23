% CONFHYP Computes the confluent hypergeometric function 
% using a series expansion:
%    f(x;a,b)= sum i (G(a)G(b+i)x^i)/(G(a+i)G(b)i),
% where G is the Gamma function.  Notice that f(x;a,a)=exp(x).
% This function solves Kummer's Equation:
%     xf''(x)+(a-x)f'(x)-bf(x)=0,
% with f(0)=0.
% The parameter n should be a scalar that determines how many terms are
% used in the series expansion.

% Copyright (c) 1997-2000, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function z=confhyp(x,a,b,n)
if nargin<4; n=100; end
z=1;
a=a+n-1; 
b=b+n-1;
while n>0
  z=z.*((a./b).*(x./n))+1;
  a=a-1; b=b-1; n=n-1;
end
return

