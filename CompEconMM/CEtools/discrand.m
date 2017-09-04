%% DISCRAND
%
%  Generates pseudo-random discrete variates.
%
%  Usage
%    j = discrand(n,p)
%  Input
%    n = number of occurances simulated
%    p = m.1 vector of probabilities with which states 1,2,...,m occur
%  Output
%    j = n.1 vector of random occurances of states 1,2,...,m

%  Copyright(c) 1997-2015
%    Mario J. Miranda - miranda.4@osu.edu
%    Paul L. Fackler  - paul_fackler@ncsu.edu

function j = discrand(n,p)
if n<1
  j = [];
else
  p = [0;cumsum(p(:))];
  m = length(p);
  m = m-length(find(p==p(end)));
  [~,j] = sort([p(1:m); rand(n,1)]);
  temp = find(j>m);
  i = j(temp)-m;
  j = temp-(1:n)';
  j(i) = j(:);
  j(j==0) = length(find(p==p(1)));
  j = max(j,1);
  j = min(j,m);
end