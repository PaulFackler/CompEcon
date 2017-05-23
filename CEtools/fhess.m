% FHESS Computes finite difference Hessian
% Uses a more flexible input/output arrangement than fdhess
% USAGE
%   H = fhess(f,ind,x1,x2,...);
% INPUTS
%   f         : name of function of form fval = f(x1,x2,...)
%   ind       : index of input/output arguments for Hessian
%   x1,x2,... : input arguments for f
% OUTPUT
%   H         : finite difference Hessian
% 
% Example: 
%   Suppose f has calling syntax: [y1,y2]=f(x1,x2,x3);
%     H=fhess(f,[3,2],x1,x2,x3)
%   computes an approximation for d^2y2/dx3^2.
% Let m=prod(size(y2)) and n=prod(size(y3)).
% The size of the output will be m by n by n, 
%   unless m=1, when the output is n by n.
%
% USER OPTIONS (SET WITH OPSET)
%   tol      : a factor used in setting the step size
%               increase if f is inaccurately computed
%               (default: eps^(1/4))
%   diagonly : computes just the diagonal elements of the Hessian 

% Copyright (c) 1997-2002, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function H = fhess(f,ind,varargin)
  tol      = optget(mfilename,'tol',eps.^(1/4));
  diagonly = optget(mfilename,'diagonly',0); 

  % Determine indices of input and output variables
  switch length(ind)
  case 2
    in=ind(1); out=ind(2);
  case 1
    in=ind; out=1;
  case 0
    in=1;  out=1;
  otherwise
    error('Incorrectly specified index variable (ind)')
  end

  % Check function inputs and outputs (if f is an m-file)
  i=exist(f);
  if i==2 | i==6         % input is an m-file or p-file
    if nargin(f)<in
      error('Index for input argument exceeds number of function inputs')
    end
    if nargout(f)<out
      error('Index for output argument exceeds number of function outputs')
    end
  end
 
  fval=cell(out,1); 
  [fval{:}] = feval(f,varargin{:});
  fx=fval{out}(:);
  m=size(fx,1);
 
  x=varargin{in}(:);
  n = size(x,1);

  % Compute the stepsize (h)
  h = tol.*max(abs(x),1);
  xh = x+h;
  h = xh-x;    
  ee = sparse(1:n,1:n,h,n,n);
 
  hh=h*h';
  H=zeros(m,n,n);
 
  % Compute forward and backward steps
   f1 = zeros(m,n);
  f0 = zeros(m,n);
  for i=1:n
    varargin{in}(:) = x+ee(:,i); [fval{:}]=feval(f,varargin{:}); f1(:,i)=fval{out}(:);
    varargin{in}(:) = x-ee(:,i); [fval{:}]=feval(f,varargin{:}); f0(:,i)=fval{out}(:);
    H(:,i,i) = (f1(:,i)+f0(:,i)-2*fx)./hh(i,i);
  end
  
  % Compute double steps
  if ~diagonly 
    for i=1:n
      for j=1:n
        if i~=j
          varargin{in}(:) = x+ee(:,i)-ee(:,j); 
          [fval{:}]=feval(f,varargin{:}); fxx=fval{out}(:);
          H(:,i,j) = (f1(:,i)+f0(:,j)-fx-fxx)./hh(i,j);
        end
      end
    end
    H=(H+permute(H,[1 3 2]))./2;       % transpose the 2 & 3 dimensions
  end

  if m==1           % if output is 1-D, return an nxn matrix
    H=squeeze(H);
  end