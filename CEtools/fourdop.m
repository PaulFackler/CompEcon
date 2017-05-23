% FOURDOP Computes differential operator for Fourier functions
% USAGE
%  [D,n,a,b,parms]=fourdop(m,rho,beta,order);
%
% See also: FOURBAS, FUNDEF

% Copyright (c) 1997-2000, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function [D,n,a,b,parms]=fourdop(m,rho,beta,order)

n=m+m+1;
D=speye(n);
v=[ones(1,m);-ones(1,m)];
u=flipud(reshape(2:n,2,m));
for i=order:-1:1
  D=sparse(1:n,[1;u(:)],[beta;v(:)],n,n)*D;
  if beta~=0, beta=beta-1; end
end
for i=order:-1
  D=sparse(1:n,1:n,[1/(beta+1);-v(:)],n,n)*D;
  %D(:,1)=D(:,1)-d(1)/2;         % adjustment to make value at original endpoint equal 0
  beta=beta+1;
end

n=m+m+1;
a=-inf;
b=inf;
parms={m,rho,beta};
