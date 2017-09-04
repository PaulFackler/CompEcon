% MINMAX
%
%  Minimax transformation for solving NCP as rootfinding problem
%
%  Usage
%    [fnew,Jnew] = minmax(x,a,b,f,J)
%  Input
%    x         : n.1 evaluation point
%    a         : n.1 lower bound
%    b         : n.1 upper bound
%    f         : n.1 function value
%    J         : n.n Jacobian
%  Output
%    fnew      : n.1 transformed function value
%    Jnew      : n.n transformed Jacobian

%  Copyright(c) 1997-2015
%   Mario J. Miranda - miranda.4@osu.edu
%   Paul L. Fackler  - paul_fackler@ncsu.edu

function [fnew,Jnew] = minmax(x,a,b,f,J)

if length(a)==1, a = a*ones(size(x)); end
if length(b)==1, b = b*ones(size(x)); end

da = a-x;
db = b-x;
fnew = min(max(f,da),db); 

if nargout==2
   Jnew = -eye(length(x));
   i = find(f>da & f<db);
   Jnew(i,:) = J(i,:);
end