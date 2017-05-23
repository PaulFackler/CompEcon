% TRIGAMMA Calculates the value of the trigamma function
% The trigamma function is the second derivative of the loggamma function.
% USAGE
%    p=trigamma(x);
% Accepts positive matrices.

% Copyright (c) 1997-2000, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function p=trigamma(x)
  x=x+6;
  p=1./(x.*x);
  p=(((((0.075757575757576*p-0.033333333333333).*p+0.0238095238095238) ...
     .*p-0.033333333333333).*p+0.166666666666667).*p+1)./x+0.5*p;
  for i=1:6
    x=x-1;
    p=p+1./(x.*x);
  end
  p(x<=0)=NaN;