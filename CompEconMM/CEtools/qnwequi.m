%% QNWEQUI
%
%  Generates equidistributed nodes and equal weights on hypercube
%
%  Usage
%    [x,w] = qnwsimp(n,a,b,type)
%  Let
%    d = dimension of integration interval
%  Input
%    n   : number of nodes and weights
%    a   : 1.d lower bounds of integration interval
%    b   : 1.d upper bounds of integration interval
%    type: type of sequence
%           N - Neiderreiter (default)
%           W - Weyl
%           H - Haber
%           R - pseudo Random
%  Output
%    x   : n.d equidistriuted nodes
%    w   : n.1 equal weights
%  Note
%    To compute definte integral of a real-valued function f defined on a
%    hypercube [a,b] in R^d, write a Matlab function f that returns an m.1 
%    vector when passed an m.n matrix, and write [x,w]=qnwequi(n,a,b,type); 
%    intf=w'*f(x).

%  Copyright(c) 1997-2015
%   Mario J. Miranda - miranda.4@osu.edu
%   Paul L. Fackler  - paul_fackler@ncsu.edu

function [x,w] = qnwequi(n,a,b,type)

if any(a>b)
  error('In qnwequi: upper bounds must exceed lower bounds.') 
end
if any(n<2)
  error('In qnwequi: number of nodes n must exceed one.')
end

global equidist_pp

if isempty(equidist_pp)
  equidist_pp = sqrt(primes(7920));   % good for d<=1000
end

d = max(length(n),max(length(a),length(b)));
n = prod(n);
if nargin<4
  type = 'N';
end

i = (1:n)';
switch upper(type(1))
  case 'N'                 % Neiderreiter
    j = 2.^((1:d)/(d+1));
    x = i*j;
    x = x-fix(x);
  case 'W'                 % Weyl
    j = equidist_pp(1:d);
    x = i*j;
    x = x-fix(x);
  case 'H'                 % Haber
    j = equidist_pp(1:d);
    x = (i.*(i+1)./2)*j;
    x = x-fix(x);
  case 'R'                 % pseudo-random
    x = rand(n,d);
  otherwise
    error('In qnwequi: unknown sequence requested.')
end

u = ones(n,1);
r = b-a;
x = a(u,:) + x.*r(u,:);
w = (prod(r)/n)*u;