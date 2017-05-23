% RICCATI Solves r0+r1C+r11Cr12+r21Cr22C=0
% This Riccati equation arises in RE models
% USAGE
%   [c,e]=riccati(c,r0,r1,r11,r12,r21,r22)
% INPUTS
%   c                     : m.n matrix of starting values
%   r0,r1,r11,r12,r21,r22 : parameter values (m.n, m.m, m.m, n.n, m.n, n.m)
% OUTPUTS
%   c : solution value
% 
% Used by RELIN

% Copyright (c) 1997-2002,  Paul L. Fackler & Mario J. Miranda 
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function [e,J,c]=riccati(c,r0,r1,r11,r12,r21,r22)
maxit=100;
tol=1e-12;

[m,n]=size(r0);
R0=r0(:);
R1=kron(eye(n),r1)+kron(r12',r11);

if isempty(c)
 c=R1\(-R0);
else
  c=c(:);
end

for i=1:maxit
  temp1=r21*reshape(c,m,n);
  temp2=r22*reshape(c,m,n);
  temp=temp1*temp2;
  e=R0+R1*c+temp(:);
  if max(abs(e))<tol, break; end
  temp=kron(temp2',r21);
  J=(R1+temp+kron(eye(n),temp1*r22));
warning off
  dc=J\e;
warning on
  c=c-dc;
end

if i>=maxit,
  disp(['Procedure did not converge in ' num2str(i) ' iterations'])
end

c=reshape(c,m,n);
e=reshape(e,m,n);
