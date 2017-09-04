%% FDJAC
%
%  Computes two-sided finite difference Jacobian of function from R^n to R^m
%
%  Usage
%    J = fdjac(f,x,varargin)
%  Input
%    f         : R^m-valued function of form fval=f(x,varargin) where fval 
%                is analytically computed value of f
%    x         : n.1 evaluation point
%    varargin  : optional parameters passed to f
%  Output
%    J         : m.n finite difference Jacobian
%  Options
%    tol       : factor used to set step size (eps^(1/3))

%  Copyright(c) 1997-2015
%   Mario J. Miranda - miranda.4@osu.edu
%   Paul L. Fackler  - paul_fackler@ncsu.edu

function J = fdjac(f,x,varargin)

% Set options to defaults, if not set by user
tol = optget(mfilename,'tol',eps.^(1/3));

% Compute stepsize
h = tol.*max(abs(x),1);
xh1 = x+h;
xh0 = x-h;
h   = xh1-xh0;

% Compute finite difference Jacobian
for j=1:length(x);
  xx = x;
  xx(j) = xh1(j); f1=feval(f,xx,varargin{:});
  xx(j) = xh0(j); f0=feval(f,xx,varargin{:});
  J(:,j) = (f1-f0)/h(j);
end

% One-sided derivative
% % Compute stepsize 
% h = tol.*max(abs(x),1);
% xh1 = x+h; 
% xh0 = x;
% h   = xh1-xh0;
% 
% % Compute finite difference Jacobian
% for j=1:length(x);
%   xx = x;
%   xx(j) = xh1(j); f1=feval(f,xx,varargin{:});
%   xx(j) = xh0(j); f0=feval(f,xx,varargin{:});
%   J(:,j) = (f1-f0)/h(j);
% end