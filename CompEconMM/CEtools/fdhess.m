%% FDHESS 
%
%  Computes two-sided finite difference Hessian of real-valued function on R^n
%
%  Usage
%    H = fdhess(f,x,varargin)
%  Input
%    f         : real-valued function of form fval=f(x,varargin) where fval 
%                is analytically computed value of f
%    x         : n.1 evaluation point
%    varargin  : optional parameters passed to f
%  Output
%    H         : n.n finite difference Hessian
%  Options
%    tol       : factor used to set step size (eps^(1/4))
%    diagonly  : compute only diagonal elements of Hessian (0)

%  Copyright(c) 1997-2015
%   Mario J. Miranda - miranda.4@osu.edu
%   Paul L. Fackler  - paul_fackler@ncsu.edu

function H = fdhess(f,x,varargin)

% Set options to defaults, if not set by user
tol      = optget(mfilename,'tol',eps.^(1/4));
diagonly = optget(mfilename,'diagonly',0);

k  = size(x,1);
fx = feval(f,x,varargin{:});

% Compute stepsize
h  = tol*max(abs(x),1);
xh = x+h;
h  = xh-x;
ee = sparse(1:k,1:k,h,k,k);

% Compute forward and backward steps
gplus  = zeros(k,1);
gminus = zeros(k,1);
for i=1:k
  gplus(i)  = feval(f,x+ee(:,i),varargin{:});
  gminus(i) = feval(f,x-ee(:,i),varargin{:});
end

H=h*h';
% Compute double steps
if diagonly
  for i=1:k
    H(i,i) = (gplus(i)+gminus(i)-2*fx)/ H(i,i);
  end
else
  for i=1:k
    for j=1:k
      if i==j
        H(i,j) = (gplus(i)+gminus(j)-2*fx)/ H(i,j);
      else
        fxx=feval(f,x+ee(:,i)-ee(:,j),varargin{:});
        H(i,j) = (gplus(i)+gminus(j)-fx-fxx)/ H(i,j);
      end
    end
  end
  H=(H+H')/2;
end