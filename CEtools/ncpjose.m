% NCPJOSE  Solves nonlinear complementarity problem using sequential LCP method
%     a <= x <= b
%     x_i > a_i => f_i(x) => 0
%     x_i < b_i => f_i(x) =< 0
% USAGE
%   [x,fval] = ncpjose(lcp,f,a,b,x,varargin);
% INPUTS
%   lcp     : name of linear complementarity solver
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
%   tol     : convergence tolerance
%   maxit   : maximum number of iterations

% Copyright (c) 1997-2000, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function [x,fval] = ncpjose(lcp,f,a,b,x,varargin)

maxit    = optget('ncpjose','maxit',100);
tol      = optget('ncpjose','tol',sqrt(eps));


optset('lcpsolve','maxit',20);
warning off

for it=1:maxit
   xold = x;
   [fval,fjac] = feval(f,x,varargin{:});
   if norm(minmax(x,a,b,fval),inf)<tol, return; end
   switch lcp
   case 'smooth'
      optset('lcpsolve','type','smooth');
      [x,z] = lcpsolve(fjac,fval-fjac*xold,a,b,xold);
   case 'minmax'
      optset('lcpsolve','type','minmax');
      [x,z] = lcpsolve(fjac,fval-fjac*xold,a,b,xold);
   case 'baard'
      [x,z] = lcpbaard(fjac,fval-fjac*xold,a,b,xold);
   case 'lemke'
      [x,z] = lcplemke(fjac,fval-fjac*xold,a,b,xold);
   end
   if any(isnan(x)) || any(isinf(x))
     warning('NaNs or Infs encountered - terminating iterations')
     break
   end
end

optset('lcpsolve','defaults');
warning on

warning('Failure to converge in ncpjose');
fval = feval(f,x,varargin{:});

