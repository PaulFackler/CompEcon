%% CHECKJAC 
%
%  Compares analytic and finite difference Jacobian of function from R^n to R^m
%
%  Usage
%    [error,i,j] = checkjac(f,x,varargin)
%  Input
%    f         : function of form [fval,fjac]=f(x) where fval and fjac are 
%                analytically computed value and Jacobian of f
%    x         : evaluation point
%    varargin  : optional parameters passed to f
%  Output
%    error     : maximum discrepancy between analytic and finite difference Jacobians
%    i,j       : indices of maximum discrepancy

%  Copyright(c) 1997-2015
%   Mario J. Miranda - miranda.4@osu.edu
%   Paul L. Fackler  - paul_fackler@ncsu.edu

function [error,i,j] = checkjac(f,x,varargin)

[~,fjac]  = feval(f,x,varargin{:});
if isempty(fjac)
  return
end
fjacfd    = fdjac(f,x,varargin{:});
[error,i] = max(abs(fjac-fjacfd));
[error,j] = max(error);
i = i(j);
if error>1.e-6
  fprintf('WARNING: Analytic Jacobian may be incorrectly coded\n')
  fprintf('  Maximum Discrepancy %12.6e\n',error)
  fprintf('  Row     %7i\n',i)
  fprintf('  Column  %7i\n',j)
end