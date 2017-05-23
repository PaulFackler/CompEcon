% FUNDOP Computes derivative operators
% USAGE 
%  [D,Dfspace]=fundop(fspace,order);
% 
% To create a function (family definition structure, Dfspace, and coefficient matrix, Dc)
% equal to the derivative of integral of an existing function (fspace and c) use:
%   [D,Dfspace]=fundop(fspace,order);
%   Dc=D*c;

% Copyright (c) 1997-2002, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function [D,Dfspace]=fundop(fspace,order);

if length(order)==1, order=order(ones(fspace.d,1)); end

D=cell(1,fspace.d);
Dfspace=fspace;
for i=1:fspace.d
  if order(i)~=0
    [d,n,a,b,parms]=feval([fspace.bastype{i} 'dop'],fspace.parms{i}{:},order(i));
    D{i}=d{end};
    Dfspace.n(i)=n;
    Dfspace.a(i)=a;
    Dfspace.b(i)=b;
    Dfspace.parms{i}=parms;
  else
    D{i}=speye(fspace.n(i));
  end
end

if fspace.d>1
  D=ckron(D(fspace.d:-1:1));
else
  D=D{1};
end