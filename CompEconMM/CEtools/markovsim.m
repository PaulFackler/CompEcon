%% MARKOVSIM
%
%  Simulates transitions of a Markov chain.
%
%  Usage
%    j = markovsim(i,Q,p)
%  Input
%    i   n.1  integers between 1 & m
%    Q   m.m  probability transition matrix
%    p   optional number of periods to be simulated (default is 1)
%  Output
%    j = n.p  matrix of integers between 1 & m, generated according to Q

%  Copyright(c) 1997-2015
%   Mario J. Miranda - miranda.4@osu.edu
%   Paul L. Fackler  - paul_fackler@ncsu.edu

function j = markovsim(i,Q,p)
m = size(Q,1);
n = length(i);
if nargin<3
  r = rand(n,1);
  j = min(sum(r(:,ones(1,m))>cumsum(Q(i,:),2),2)+1,m);
else
  j = zeros(n,p);
  for k=1:p
      r = rand(n,1);
      i = min(sum(r(:,ones(1,m))>cumsum(Q(i,:),2),2)+1,m);
      j(:,k) = i;
  end
end