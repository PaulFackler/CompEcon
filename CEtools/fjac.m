% FJAC Computes two-sided finite difference Jacobian
% Uses a more flexible input/output arrangement than fdjac
% USAGE
%   J = fjac(f,ind,x1,x2,...)
% INPUTS
%   f         : name of function of form fval = f(x1,x2,...)
%   ind       : index of input/output arguments for Jacobian
%   x1,x2,... : input arguments for f
% OUTPUT
%   J         : finite difference Jacobian
%
% Example: 
%   Suppose f has calling syntax: [y1,y2]=f(x1,x2,x3);
%     J=fjac(f,[3,2],x1,x2,x3)
%   computes an approximation for dy2/dx3.
% Let m=prod(size(y2)) and n=prod(size(y3)).
% The size of the output will be m by n. 
%
% USER OPTIONS (SET WITH OPSET)
%   tol      : a factor used in setting the step size
%               change if f is inaccurately computed 
%               (default: eps^(1/3))

% Copyright (c) 1997-2002, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function J = fjac(f,ind,varargin)
tol = optget(mfilename,'tol',eps.^(1/3));

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
if ~isa(f,'function_handle')
i=exist(f);
if i==2 | i==6         % input is an m-file or p-file
  if nargin(f)<in & nargin(f)>=0
    error('Index for input argument exceeds number of function inputs')
  end
  if nargout(f)<out & nargout(f)>=0
    error('Index for output argument exceeds number of function outputs')
  end
end
end

x=varargin{in}(:);
h = tol.*max(abs(x),1);
xh1=x+h; xh0=x-h;
h=xh1-xh0;
fval=cell(out,1);
for j=1:length(x);
  varargin{in}(:) = x; 
  varargin{in}(j) = xh1(j); [fval{:}]=feval(f,varargin{:}); f1=fval{out};
  varargin{in}(j) = xh0(j); [fval{:}]=feval(f,varargin{:}); f0=fval{out};
  % preallocate output when j=1
  if j==1, J=zeros(size(f1(:),1),size(x(:),1)); end  
  J(:,j) = (f1(:)-f0(:))/h(j);
end