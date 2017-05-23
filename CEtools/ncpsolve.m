% NCPSOLVE  Solves nonlinear complementarity problem
%       a <= x <= b
%     x_i > a_i => f_i(x) => 0
%     x_i < b_i => f_i(x) =< 0
% Uses either minimax or semismooth safeguarded Newton method
% USAGE
%   [x,fval] = ncpsolve(f,a,b,x,varargin)
% INPUTS
%   f       : name of function of form [fval,fjac]=f(x)
%   a       : n by 1 vector, left bound on x
%   b       : n by 1 vector, right bound on x
%   x       : n by 1 vector, initial guess for solution
%   varargin: additional arguments for f (optional)
% OUTPUTS
%   x       : solution to ncp
%   fval    : function value at x
%
% Setable options (use OPTSET):
%   tol       : convergence tolerance
%   maxit     : maximum number of iterations
%   maxsteps  : maximum number of backsteps
%   showiters : display results of each iteration
%   type      : rootproblem transform, 'smooth' or 'minmax'

% Copyright (c) 1997-2000, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function [x,fval] = ncpsolve(f,a,b,x,varargin)

maxit     = optget('ncpsolve','maxit',100);
tol       = optget('ncpsolve','tol',sqrt(eps));
maxsteps  = optget('ncpsolve','maxsteps',10);
showiters = optget('ncpsolve','showiters',0);
type      = optget('ncpsolve','type','smooth');

if nargin < 4, x=zeros(size(a)); end

for it=1:maxit
   [fval,fjac] = feval(f,x,varargin{:});
   [ftmp,fjac] = feval(type,x,a,b,fval,fjac);
   fnorm = norm(ftmp,inf);
   if fnorm<tol, return, end
   dx = -(fjac\ftmp);
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
   if showiters, fprintf('%4i %4i %6.2e\n',[it backstep fnormnew]); end
end
warning('Failure to converge in ncpsolve')
fval = feval(f,x,varargin{:});
