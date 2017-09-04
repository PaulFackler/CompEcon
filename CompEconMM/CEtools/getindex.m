%% GETINDEX 
%
%  Finds row of a matrix that most reseemling a vector or series of vectors.
%
%  Usage
%    i = getindex(s,S);
%  Input
%    s         : 1.d vector or p.d matrix
%    S         : n.d matrix
%  Output
%    i         : p.1 vector of integers in {1,...,n} indicating row of S
%                that most closely matches each row in s
%  Options
%  Example
%     S=[0 0; 0 1; 1 0; 1 1]; s=[0 0; 0 0; 1 0; 0 0; 1 1]; getindex(s,S)
%     returns [1; 1; 3; 1; 4]

%  Copyright(c) 1997-2015
%   Mario J. Miranda - miranda.4@osu.edu
%   Paul L. Fackler  - paul_fackler@ncsu.edu

function i = getindex(s,S)

[n,d] = size(S);
p = size(s,1);
if p==1
  [junk,i] = min(sum(abs(s(ones(n,1),:)-S),2));
else
  [S,s] = gridmake(S,s);
  z = sum(abs(s-S),2);
  [junk,i] = min(reshape(z,n,p));
  i = i';
end