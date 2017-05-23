% AFFASSET Solves affine asset pricing models
%  dX/dt = A'X + 0.5B'diag(C'X)C'X-g
%  dx/dt = a'X + 0.5b'diag(C'X)C'X-g0
% USAGE
%   [X,x]=affasset(t,a,A,b,B,C,g,G,h,h0);
% INPUTS
%   t              : a m-vector of time values
%   a,A,b,B,C,g,g0 : parameters of n-state affine model
%   h,h0           : initial values h=X(t(1)) and h0=x(t(1)) 
% OUTPUTS
%   X : m by n solution matrix for X(t)
%   x : m by 1 solution vector for x(t)
%
% Uses: affode
%
% Note: CompEcon text refers to X and x as beta and beta0, respectively.

% Copyright (c) 1997-2002, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function [X,x]=affasset(t,a,A,b,B,C,g,g0,h,h0);
% Define parameters for ODE file (merges equations for X and x)
AA=[A.';a.'];
BB=[B.';b.']/2;
GG=[g(:);g0];

%affode = @(t,x)

% Call solver
[t,X]=ode45('affode',t,[h(:);h0],[],AA,BB,C.',GG);

% Break apart results
x=X(:,end);
X(:,end)=[];
