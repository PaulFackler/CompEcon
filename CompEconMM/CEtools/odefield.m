%% ODEFIELD
%
%  Generates velocity field for 2-dimensional ODE
%
%  Specifically, generates velocity field for 2-dimensionnal 1st-order ODE
%    x'(t) = f(x(t)), t in [0,T]
%  Here, x is 2.1 vector-valued function defined on time domain [0,T] and
%  x' is its 2.1 vector-valued derivative with respect to time.
%
%  Usage
%   odefield(f,x1lim,x2lim,n,varargin)
%  Input
%    f         : velocity function (see below)
%    x1lim     : 2.1 lower & upper limits of x1
%    x2lim     : 2.1 lower & upper limits of x2
%    n         : 2.1 number of coordinates per dimension
%    varargin  : optional parameters passed to f
%  Output
%    Generates n(1) by n(2) velocity field, assuming figure has
%    already been defined by calling program.

%  Copyright(c) 1997-2015
%   Mario J. Miranda - miranda.4@osu.edu
%   Paul L. Fackler  - paul_fackler@ncsu.edu

function odefield(f,x1lim,x2lim,n,varargin)
if nargin<5, n=[]; end
if isempty(n)
  n = [11 11];
end
x = gridmake(nodeunif(n(1),x1lim(1),x1lim(2)),nodeunif(n(2),x2lim(1),x2lim(2)))';
v = real(f(x,varargin{:}));
quiver(x(1,:),x(2,:),v(1,:),v(2,:),'o','LineWidth',1,'MarkerSize',3,'Color',[0.6 0.6 0.6]);