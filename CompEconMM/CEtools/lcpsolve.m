%% LCPSOLVE 
%
%  Solves linear complementarity problem using Newton's method
%
%  Problem
%     z = M*x + q
%     a <= x <= b
%     x_i > a_i => z_i => 0
%     x_i < b_i => z_i =< 0
%
%  Usage
%    [x,z] = lcpsolve(M,q,a,b,x)
%  Input
%    M         : n.n matrix
%    q         : n.1 vector
%    a         : n.1 lower bound on x
%    b         : n.1 upper bound on x
%    x         : n.1 initial guess for solution
%  Output
%    x         : solution to lcp
%    z         : function value at x
%  Options
%    maxit     : maximum number of iterations (100)
%    tol       : convergence tolerance (sqrt(eps))
%    maxsteps  : maximum number of backsteps (20)
%    type      : rootproblem transform, 'ssmooth' or 'minmax' (ssmooth)
%    showiters : display results of each iteration (1)

%  Copyright(c) 1997-2015
%   Mario J. Miranda - miranda.4@osu.edu
%   Paul L. Fackler  - paul_fackler@ncsu.edu

function [x,z] = lcpsolve(M,q,a,b,x)

maxit     = optget('lcpsolve','maxit',100);
tol       = optget('lcpsolve','tol',sqrt(eps));
maxsteps  = optget('lcpsolve','maxsteps',20);
type      = optget('lcpsolve','type','ssmooth');
showiters = optget('lcpsolve','showiters',1);

if showiters, fprintf('\nIn lcpsolve:\n'), end

if nargin < 5, x=q; end

for it=1:maxit
  z = M*x+q;
  [ztmp,Mtil] = feval(type,x,a,b,z,M);
  fnorm = norm(ztmp);
  if showiters, fprintf('%4i %6.2e\n',[it fnorm]); end
  if fnorm<tol, return, end
  dx = -(Mtil\ztmp);
  fnormold = inf;
  for backstep=1:maxsteps
    xnew = x + dx;
    znew = M*xnew+q;
    ztmp = feval(type,xnew,a,b,znew);
    fnormnew = norm(ztmp);
    if fnormnew<fnorm, break, end
    if fnormold<fnormnew, dx=dx*2; break, end
    fnormold = fnormnew;
    dx = dx/2;
  end
  x = x+dx;  
end
warning('In lcpsolve" failure to converge.')
z = M*x+q;