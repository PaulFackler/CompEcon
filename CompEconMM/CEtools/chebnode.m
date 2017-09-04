%% CHEBNODE
%
%  Computes standard nodes for Chebychev polynomials
%
%  Usage
%    x = chebnode(n,a,b)
%  Input
%    n         : number of nodes
%    a         : lower bound of approximation interval
%    b         : upper bound of approximation interval
%  Output
%    x         : n.1 vector of nodes
%  Options
%    nodetype  0 : Gaussian nodes (do not include endpoints)
%              1 : Gaussian nodes extended to endpoints
%              2 : Lobatto nodes  (includes endpoints)

%  Copyright(c) 1997-2015
%   Paul L. Fackler  - paul_fackler@ncsu.edu
%   Mario J. Miranda - miranda.4@osu.edu


function x = chebnode(n,a,b)

nodetype = optget('chebnode','nodetype',0);

s = (b-a)/2;
m = (b+a)/2;

if nodetype<2                           % usual nodes
  k = pi*(0.5:(max(n)-0.5))';
  x = m(1)-cos(k(1:n(1))/n(1))*s(1);
  if nodetype==1                        % Extend nodes to endpoints
    aa = x(1);
    bb = x(end);
    x = (bb*a-aa*b)/(bb-aa)+(b-a)/(bb-aa)*x;
  end
else                                    % Lobatto nodes
  k = pi*(0:(max(n)-1))';
  x = m(1)-cos(k(1:n(1))/(n(1)-1))*s(1);
end