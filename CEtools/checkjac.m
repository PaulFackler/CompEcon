% CHECKJAC Compares analytic and finite difference derivative
% USAGE
%   [error,i,j] = checkjac(f,x,varargin);
% INPUTS
%   f       : name of function of form [fval,fjac]=f(x)
%             where fval and fjac are value and jacobian of f
%   x       : evaluation point
%   varargin: additional arguments for f (optional)
% OUTPUTS
%   error   : maximum difference between analytical 
%               and finite difference Jacobians
%   i,j     : indices of maximum difference
%
% USES: fdjac

% Copyright (c) 1997-2000, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function [error,i,j] = checkjac(f,x,varargin)

[fval,fjac] = feval(f,x,varargin{:});
fjacfindif = fdjac(f,x,varargin{:});
[error,i] = max(abs(fjac-fjacfindif));
[error,j] = max(error);
i = i(j);
