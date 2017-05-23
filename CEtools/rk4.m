% RK4 Solves initial value problems using fourth-order Runge-Kutta.
%       X' = F(X) with X(A)=X0.
% where X is a D-dimensional vector.
% USAGE
%   [T,X]=rk4(f,T,X0,options,varargin)
% INPUTS
%   F: the defining differential equation (supplied by the user).
%      F must take three inputs:
%         T: the current time value
%         X: where X is DxK
%         Flag: not used - for compatability with MATLAB's solvers
%         Optional additional parameters
%      F returns a DxK matrix.
%   T: a vector of N time values
%   X0: a DxK matrix of initial values
%   options: not used -  for compatability with MATLAB's solvers
%   Optional additional parameters to pass to F
% OUTPUTS
%   T: Same as input value - for compatability with MATLAB's solvers
%   X: An NxDxK matrix with dimensions
%      1: time
%      2: variable
%      3: starting point
%
% Comments:  this program is written so the system can be solved at K initial
% values simultaneously.  The user supplied function F should therefore be
% written to accept matrix values of X.  This feature makes it possible
% to solve for an entire phase diagram with one call to RK4.
% The vector of time values should be sorted (ascending or descending).
% The values need not be equally spaced; this allows smaller steps to 
% be used during times of non-smooth behavior (assuming one knows 
% when this occurs a priori).

% Copyright (c) 1997-2000, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function [T,X]=rk4(f,T,X0,options,varargin)

  [d,k]=size(X0);
  n=length(T);
  X=zeros(n,d,k);                     % initialize X
  X(1,:,:)=X0;
  h=[0;diff(T)];                      % Determine size of time steps
  t=T(1);
  for i=2:n
    hh=h(i);
    f1=feval(f,t,X0,'',varargin{:})*(hh/2);
    f2=feval(f,t+hh/2,X0+f1,'',varargin{:})*hh;
    f3=feval(f,t+hh/2,X0+f2/2,'',varargin{:})*hh;
    f4=feval(f,t+hh,X0+f3,'',varargin{:})*(hh/2);
    X0=X0+(f1+f2+f3+f4)/3;
    t=t+hh;
    X(i,:,:)=X0;
  end
