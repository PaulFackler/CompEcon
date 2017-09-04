%% GOLDEN
%
%  Computes local maximum of univariate function on interval using Golden Search
%
%  Usage
%    [x,fval] = golden(f,a,b,varargin)
%  Input
%    f         : real-valued function of form fval=f(x)
%    a         : lower bound of interval
%    b         : upper bound of interval
%    varargin  : optional parameters passed to f
%  Output
%    x         : local maximum of f
%    fval      : function value at x
%  Options
%    tol       : convergence tolerance

%  Copyright(c) 1997-2015
%   Mario J. Miranda - miranda.4@osu.edu
%   Paul L. Fackler  - paul_fackler@ncsu.edu

function [x1,f1] = golden(f,a,b,varargin)

tol = optget('golden','tol',sqrt(eps));

alpha1 = (3-sqrt(5))/2;
alpha2 = (sqrt(5)-1)/2;
d  = b-a;
x1 = a+alpha1*d;
x2 = a+alpha2*d;
f1 = feval(f,x1,varargin{:});
f2 = feval(f,x2,varargin{:});

d = alpha1*alpha2*d;
while d>tol
  d = d*alpha2;
  if f2<f1 % x2 is new upper bound
    x2 = x1; x1 = x1-d;
    f2 = f1; f1 = feval(f,x1,varargin{:});
  else     % x1 is new lower bound
    x1 = x2; x2 = x2+d;
    f1 = f2; f2 = feval(f,x2,varargin{:});
  end
end

% Return larger of the two
if f2>f1, x1 = x2; f1 = f2; end