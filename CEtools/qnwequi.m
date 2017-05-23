% QNWEQUI Generates equidistributed sequences
%         with property that averages value of
%         integrable function evaluated over the sequence
%         converges to the integral as n goes to infinity
% USAGE
%   [x,w] = qnwequi(n,a,b,type);
% INPUTS
%   n   : the number of sequence points
%   a,b : d-vectors of left and right endpoints for each dimension
%        (problem dimension, d, is determined by a and b)
%   type: type of sequence
%         N - Neiderreiter (default)
%         W - Weyl
%         H - Haber
%         R - pseudo Random
% OUTPUTS
%   x   : an nxd matrix of equidistributed (or uniformly random) nodes
%   w   : n by 1 vector of weights
%
% Note: for d>1000 using type='W' or type='H' modify the
%   definition of pp to use a larger prime sequence.
%
% Reference: Judd (sec. 9.1) for discussion

% Copyright (c) 1997-2000, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function [x,w] = qnwequi(n,a,b,type)

global equidist_pp

if isempty(equidist_pp)
  equidist_pp=sqrt(primes(7920));   % good for d<=1000 
end

d  = max(length(n),max(length(a),length(b)));
n=prod(n);
if nargin<4, type='N'; end

i=(1:n)';
switch upper(type(1))
  case 'N'                 % Neiderreiter 
    j=2.^((1:d)/(d+1));
    x=i*j;
    x=x-fix(x);
  case 'W'                 % Weyl
    j=equidist_pp(1:d);
    x=i*j;
    x=x-fix(x);
  case 'H'                 % Haber
    j=equidist_pp(1:d);
    x=(i.*(i+1)./2)*j;
    x=x-fix(x);
  case 'R'                 % pseudo-random
    x=rand(n,d);
  otherwise
    error('Unknown sequence requested')
end

u=ones(n,1);
r = b-a;
x = a(u,:) + x.*r(u,:);
w = (prod(r)/n)*u;
