% VechInv Reverses the vech operator
% USAGE:
%   A=vechinv(x,flag)
% INPUTS:
%   x : a vector of length d(d+1)/2 for some integer d
%   flag : an optional parameter (see below)
% OUTPUT:
%  A  : an nxn symmetric or triangular, depending on FLAG.
%
%  If FLAG = 0: A is a symmetrix matrix.
%            1: A is an upper triangular matrix.
%            2: A is a lower triangular matrix.

function A=vechinv(x,flag);
if nargin<2 | isempty(flag), flag=0; end

d1=length(x);
d=round((sqrt(1+8*d1)-1)/2);
if d*(d+1)/2~=d1
  error('Input is not of proper length');
end

A=zeros(d,d);
ind=find(tril(ones(d,d)));
A(ind)=x;
switch flag
case 0
  % symmetric
  A=A+triu(A',1);
case 1
  % upper triangular
  A=A';
case 2
  % lower triangular
  % already done
otherwise
  error('Invalid flag')
end