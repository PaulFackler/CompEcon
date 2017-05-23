% FDJAC Computes two-sided finite difference Jacobian
% USAGE
%   fjac = fdjac(f,x,P1,P2,...)
% INPUTS
%   f         : name of function of form fval = f(x)
%   x         : evaluation point
%   P1,P2,... : additional arguments for f (optional)
% OUTPUT
%   fjac      : finite differnce Jacobian
%
% USER OPTIONS (SET WITH OPSET)
%   tol       : a factor used in setting the step size
%               increase if f is inaccurately computed

% Copyright (c) 1997-2002, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function fjac = fdjac(f,x,varargin)

tol    = optget(mfilename,'tol',eps.^(1/3));

h = tol.*max(abs(x),1);
xh1=x+h; xh0=x-h;
h=xh1-xh0;
for j=1:length(x);
   xx = x; 
   xx(j) = xh1(j); f1=feval(f,xx,varargin{:});
   xx(j) = xh0(j); f0=feval(f,xx,varargin{:});
   fjac(:,j) = (f1-f0)/h(j);
end