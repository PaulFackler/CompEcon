% FDJAC1 Computes one-sided finite difference Jacobian
% USAGE
%   fjac = fdjac(f,x,f0,P1,P2,...)
% INPUTS
%   f         : name of function of form fval = f(x)
%   x         : evaluation point
%   f0        : f(x) if previous evaluated, otherwise pass as empty
%   P1,P2,... : additional arguments for f (optional)
% OUTPUT
%   fjac      : finite difference Jacobian
%
% USER OPTIONS (SET WITH OPSET)
%   tol       : a factor used in setting the step size
%               increase if f is inaccurately computed
%
% NOTE: input syntax is different from FDJAC

% Copyright (c) 1997-2003, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function fjac = fdjac1(f,x,f0,varargin)

tol    = optget(mfilename,'tol',sqrt(eps));

h = tol.*max(abs(x),1);
xh=x+h; 
h=xh-x;
n=length(x);

if nargin<3 | isempty(f0)
  f0=feval(f,x,varargin{:});
end

fjac=zeros(size(f0,1),n);       % preallocate memory
for j=1:n;
   xx = x; 
   xx(j) = xh(j); f1=feval(f,xx,varargin{:});
   fjac(:,j) = (f1-f0)/h(j);
end

return