% LQSOLVE Solves linear-quadratic DP models
% The LQ problem is:
%   max_x fs*s +  fx*x + 0.5s'*fss*s + x'*fxs*s +0.5x'*fxx*x
%   s.t. s' = g + gs*s + gx*x + e
% The solution is
%    x = xstar + Gamma(s-sstar)
% The shadow price (marginal value) function is
%    p = pstar + Delta(s-sstar)
% The controlled state process is
%    s' = sstar + G(s-sstar) + e
% USAGE
%   [Gamma,Delta,sstar,xstar,pstar] = lqsolve(fs,fx,fss,fxs,fxx,g,gs,gx,delta);
% INPUTS
%   fs,fx,fss,fxs,fxx : objective function parameters 
%   g,gs,gx           : state transition function parameters
%   delta             : discount factor
% OUTPUTS
%   Gamma   : response of control to state
%   Delta   : response of shadow price to state
%   G       : optimal state transition response
%
% Two algorithms are implemented, QZ and function iteration. 
% Set option 'qz' equal to 0 to use function iteration.

% Copyright (c) 1997-2002, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function [Gamma,Delta,G,sstar,xstar,pstar] = lqsolve(fs,fx,fss,fxs,fxx,g,gs,gx,delta);
qz = optget('lqsolve','qz',1);

ds = size(fs,2); 
dx = size(fx,2); 

if qz % Use QZ decomposition
  A=[eye(ds)          zeros(ds,dx+ds); 
     zeros(dx,ds+dx) -delta*gx'      ;
     zeros(ds,ds+dx)  delta*gs'      ];
  B=[ gs   gx   zeros(ds,ds);
      fxs  fxx  zeros(dx,ds);
     -fss -fxs' eye(ds)     ];
  [S,T,Q,Z]=qzordered(A,B);
  C=real(Z(ds+1:end,1:ds)/Z(1:ds,1:ds));
  Gamma=C(1:dx,:);
  Delta=C(dx+1:end,:);
  G=gs+gx*Gamma;
  ss=Z*((S-T)\(Q*[g(:);fx(:);-fs(:)]));
  sstar=ss(1:ds);
  xstar=ss(ds+1:ds+dx);
  pstar=ss(ds+dx+1:end);
else % Function iteration
  tol=1.e-8;
  % Starting values
  if cond(gx'*gx)>1.e10
    Delta = eye(ds);
  else
    Delta = -((delta*gx'*gx)...
             \(((delta*gs'*gx)*(delta*gx'*gs))+delta*gs'*gs-fxx'));
  end
  % Iterate on Ricatti equation
  for it=1:200
    Deltaold = Delta;
    temp  = delta*gx'*Delta*gs+fxs;
    Gamma = -((delta*gx'*Delta*gx+fxx)\temp);
    Delta = fss+delta*gs'*Delta*gs+temp'*Gamma;
    if norm(Delta-Deltaold)<tol, break, end;
  end
  G=gs+gx*Gamma;
  sstar=(eye(ds)-G)\g;
  temp=inv(delta*gx'*Delta*gx+fxx);
  xstar=-(eye(dx)+temp*(delta*gx'))\(temp*(fx'+delta*gx'*Delta*g));
  temp=(delta*gs'*Delta*gx+fxs')*temp;
  pstar=(eye(ds)-delta*gs'*Delta+temp*(delta*gx'*Delta*g+fx'))\(fs'+delta*gs'*Delta*g-temp*(fx'+delta*gx'*Delta*g));
  [Gamma,Delta,G,sstar xstar pstar]
end
