%% NCPSOLVE
%
%  Solves nonlinear complementarity problem using Newton's method
%
%  Find x in R^d such that
%     a_i <= x_i <= b_i, i=1,2,...,d
%     x_i > a_i => f_i(x) => 0, i=1,2,...,d
%     x_i < b_i => f_i(x) =< 0, i=1,2,...,d
%  where f is a function from R^d to R^d and a and b are d by 1 vectors 
%  with a<=b. Problem is solved by applying a safeguarded Newton method to
%  the equivalent re-forrmulated minimax or semismooth rootfinding problem.
%
%  Usage
%    [x,fval] = ncpsolve(func,a,b,x,varargin)
%  Input
%    func      : user-supplied function (see below)
%    a         : d.1 lower bound on x
%    b         : d.1 upper bound on x
%    x         : d.1 initial guess for solution
%    varargin  : optional parameters passed to func
%  Output
%    x         : d.1 solution to ncp
%    fval      : d.1 function value at x
%  Function File
%    func.m a user-supplied function that returns the d by 1 vector of
%    values and d by d Jacobian of the function f according to the format
%      [fval,fjac] = func(x,varargin)
%    function input
%      x         : d.1 vector
%      varargin  : optional parameters passed to func
%    function output
%      fval      : d.1 vector of function values
%      fjac      : d.d matris of partial derivatives
%  Options
%    maxit     : maximum number of iterations (100)
%    tol       : convergence tolerance (sqrt(eps))
%    maxsteps  : maximum number of backsteps (10)
%    showiters : display results of each iteration (1)
%    type      : rootproblem transform, 'ssmooth' or 'minmax' (ss)

%  Copyright(c) 1997-2015
%   Mario J. Miranda - miranda.4@osu.edu
%   Paul L. Fackler  - paul_fackler@ncsu.edu

function [x,ftmp] = ncpsolve(f,a,b,x,varargin)

maxit     = optget('ncpsolve','maxit',100);
tol       = optget('ncpsolve','tol',sqrt(eps));
maxsteps  = optget('ncpsolve','maxsteps',10);
showiters = optget('ncpsolve','showiters',1);
type      = optget('ncpsolve','type','ssmooth');

if showiters, fprintf('\nIn ncpsolve:\n'), end

if nargin<4, x=zeros(size(a)); end

% checkjac(f,x,varargin{:});

for it=1:maxit
  [fval,fjac] = feval(f,x,varargin{:});
  if isempty(fjac)
    fjac = fdjac(f,x,varargin{:});
  end
  [ftmp,fjac] = feval(type,x,a,b,fval,fjac);
  fnorm = norm(ftmp,inf);
  if showiters, fprintf('%4i %6.2e\n',[it fnorm]); end
  if fnorm<tol
    if showiters, fprintf('\n'), end
    break
  end
  dx = -fjac\ftmp;
  fnormold = inf;
  for backstep=1:maxsteps
    xnew = x + dx;
    fnew = feval(f,xnew,varargin{:});
    fnew = feval(type,xnew,a,b,fnew);
    fnormnew = norm(fnew,inf);
    if fnormnew<fnorm, break, end
    if fnormold<fnormnew, dx=2*dx; break, end
    fnormold = fnormnew;
    dx = dx/2;
  end
  x = x+dx;
end
if it==maxit
  warning('Failure to converge in ncpsolve')
end
x = real(x);
x = max(a,x);
x = min(b,x);