% FDHESS Computes finite difference Hessian
% USAGE
%   H = fdhess(f,x,P1,P2,...);
% INPUTS
%   f         : name of function of form fval = f(x)
%   x         : evaluation point
%   P1,P2,... : additional arguments for f (optional)
% OUTPUT
%   H         : finite difference Hessian
%
% USER OPTIONS (SET WITH OPSET)
%   tol       : a factor used in setting the step size
%               increase if f is inaccurately computed
%   diagonly  : computes just the diagonal elements of the Hessian 

% Copyright (c) 1997-2002, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function H = fdhess(f,x,varargin)
  tol      = optget(mfilename,'tol',eps.^(1/4));
  diagonly = optget(mfilename,'diagonly',0); 
 
  k = size(x,1);
  fx = feval(f,x,varargin{:});
 
  % Compute the stepsize (h)
  h = tol*max(abs(x),1);
  xh = x+h;
  h = xh-x;    
  ee = sparse(1:k,1:k,h,k,k);
 
  % Compute forward and backward steps
  gplus = zeros(k,1);
  gminus = zeros(k,1);
  for i=1:k
    gplus(i) = feval(f,x+ee(:,i),varargin{:});
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