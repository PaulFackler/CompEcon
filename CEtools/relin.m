% RELIN Linearizes and solves an RE model at (s,x)
% USAGE
%   [C,P,D,A,B]=relin(model,s,x)
% The solution is a decision rule of the form
%   x=Cs
% an equilibrium state transition function
%   s(t+1)=Ps(t)
% and an expectation function of the form
%   z=Ds
% Also returns A and B, where
%   A[s(t+1);x(t+1)]=B[s(t);x(t)]
%
% INPUTS
%   model    : see resolve for description
%   s        : steady state value of the state variable
%   x        : steady state value of the decsions variable
%              (s and x may be obtained using ress)
%   options  : algorithm options (structure array)
%
% OPTIONS FOR RELIN
%   loglin   : 1 if log-linear approximation is desired

% Copyright (c) 1997-2002,  Paul L. Fackler & Mario J. Miranda 
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function [C,P,D,A,B]=relin(model,s,x,options)

if nargin<4, options=[]; end
getopts(options,...
        'loglin', 0);

s=s(:)';
x=x(:)';
e=model.w'*model.e;

func=model.func;
params=model.params;

n=length(s);
m=length(x);

z=feval(func,'h',s,x,[],e,s,x,params{:});
fs=fjac(func,2,'f',s,x,z,[],[],[],params{:});
fx=fjac(func,3,'f',s,x,z,[],[],[],params{:});
fz=fjac(func,4,'f',s,x,z,[],[],[],params{:});
gs=fjac(func,2,'g',s,x,[],e,[],[],params{:});
gx=fjac(func,3,'g',s,x,[],e,[],[],params{:});
hs=fjac(func,2,'h',s,x,[],e,s,x,params{:});
hx=fjac(func,3,'h',s,x,[],e,s,x,params{:});
hsn=fjac(func,6,'h',s,x,[],e,s,x,params{:});
hxn=fjac(func,7,'h',s,x,[],e,s,x,params{:});

if loglin
  fs=fs*diag(s);
  fx=fx*diag(x);
  fz=fz*diag(z);
  gs=gs*diag(s);
  gx=gx*diag(x);
  hs=hs*diag(s);
  hx=hx*diag(x);
  hsn=hsn*diag(s);
  hxn=hxn*diag(x);
end

if loglin
  A=diag(s);
else
  A=eye(n);
end

A=[    A    zeros(n,m);
   -fz*hsn    -fz*hxn ];
B=[ gs        gx       ;
    fs+fz*hs  fx+fz*hx];

% Use the QZ decomposition
if 1
  [S,T,Q,Z]=qzordered(A,B);
  C=real(Z(n+1:end,1:n)/Z(1:n,1:n));
  P=real((Z(1:n,1:n)*(S(1:n,1:n)\T(1:n,1:n)))/Z(1:n,1:n));
% Use Newton's method 
% (does not guarantee the stable solution is found)
else
  r0=fs+fz*hs+fz*hsn*gs;
  r1=fx+fz*(hx+hsn*gx);
  r11=fz*hxn;
  r12=gs;
  r21=fz*hxn;
  r22=gx;
  [e,J,C]=riccati([],r0,r1,r11,r12,r21,r22);
  P=gs+gx*C;
end

D=(hsn+hxn*C)*P+hs+hx*C;