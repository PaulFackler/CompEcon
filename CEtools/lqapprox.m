% LQAPPROX Forms and solves linear-quadratic approximation of DP model
% USAGE
%   [vlq,xlq,plq,sstar,xstar,pstar,Gamma,Delta,G] = lqapprox(model,s,sstar,xstar,pstar)
% INPUTS
%   model : dynamic program model structure (see dpsolve)
%   s     : n-vector of state values at which approximation is evaluated
%   sstar : steady-state state variable value
%   xstar : steady-state action variable value
%   pstar : steady-state costate variable value
% OUTPUTS
%   vlq   : n-vector of values of approximate value function
%   xlq   : n-row matrix of values of action variable
%   plq   : n-row matrix of values of the shadow price
%   sstar : steady-state state variable value
%   xstar : steady-state action variable value
%   pstar : steady state shadow price
%   Gamma : response of control to state
%   Delta : response of shadow price to state
%   G     : optimal state transition response
% 
% USER OPTIONS (SET WITH OPSET)
%   qz      : 1 to use the qz algorithm
%             0 to iterate on the Riccati equation
%             default: 1
%   foc     : 1 to linearize the first order conditions (uses qz algorithm)
%             0 to solve the LQ approximation
%             default: 0 
%   ss      : 1 to compute the steady state 
%               (sstar, xstar and pstar are treated as initial values)
%             default: 1
% 
% The QZ algorithm is direct and more reliable than function iteration
% An alternative to LQ approximation is to linearize the first order
% conditions (Euler equations) and to solve the resulting linear dynamic system
%
% Numerical derivatives are used to compute partial derivatives with respect to s
% To avoid the use of numerical derivatives, code the model function file
% to return 
%   f,fx,fxx,fs,fxs,fss when the 'f' flag is passed 
% and
%   g,gx,gxx,gs,gxs,gss when the 'g' flag is passed
% The sizes of these outputs are:
%      f:ns.1, fx:ns.dx, fxx:ns.dx.dx, fs:ns.ds, fxs:ns.dx.ds, fss:ns.ds.ds
%      g:ns.ds, gx:ns.ds.dx, gxx:ns.ds.dx.dx, gs:ns.ds.ds, gxs:ns.ds.dx.ds, gss:ns.ds.ds.ds
% (gxs and gss are only needed if the foc option is set to 1)
%      
% See dpsolve for further discussion

% Copyright (c) 1997-2002, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function [vlq,xlq,plq,sstar,xstar,pstar,Gamma,Delta,G] = lqapprox(model,s,sstar,xstar,pstar)

qz  = optget('lqapprox','qz',1);
foc = optget('lqapprox','foc',0);
ss  = optget('lqapprox','ss',1);
tol=1.e-8;

if isfield(model,'w'), estar=model.w'*model.e; 
else,                  estar=[];
end

ds = length(sstar); 
dx = length(xstar); 
func   = model.func;
params = model.params;
delta  = model.discount;

% COMPUTE STEADY STATE VALUES 
if ss
  [sstar,xstar,pstar]=dpss(model,sstar,xstar,pstar);
end

% GET DERIVATIVE VALUES
if foc
  qz=1;                 % Riccati method not supported with foc option
  if nargout(func)<6    % only x derivatives available
    [f,fx,fxx] = feval(func,'f',sstar,xstar,estar,params{:});
    [g,gx,gxx] = feval(func,'g',sstar,xstar,estar,params{:});
    % Compute numerical derivatives w.r.t. s
    fs  = fjac(func,2,'f',sstar,xstar,estar,params{:});
    fxs = fjac(func,[2,2],'f',sstar,xstar,estar,params{:});
    fss = fhess(func,2,'f',sstar,xstar,estar,params{:});
    gs  = fjac(func,2,'g',sstar,xstar,estar,params{:});
    gxs = fjac(func,[2,2],'g',sstar,xstar,estar,params{:});
    gss = fhess(func,2,'g',sstar,xstar,estar,params{:});
  else                  % s and x derivatives available
    [f,fx,fxx,fs,fxs,fss] = feval(func,'f',sstar,xstar,estar,params{:});
    [g,gx,gxx,gs,gxs,gss] = feval(func,'g',sstar,xstar,estar,params{:});
  end
else
  if nargout(func)<6    % only x derivatives available
    [f,fx,fxx] = feval(func,'f',sstar,xstar,estar,params{:});
    [g,gx,gxx] = feval(func,'g',sstar,xstar,estar,params{:});
    % Compute numerical derivatives w.r.t. s
    fs  = fjac(func,2,'f',sstar,xstar,estar,params{:});
    fxs = fjac(func,[2,2],'f',sstar,xstar,estar,params{:});
    fss = fhess(func,2,'f',sstar,xstar,estar,params{:});
    gs  = fjac(func,2,'g',sstar,xstar,estar,params{:});
  else                  % s and x derivatives available
    [f,fx,fxx,fs,fxs,fss] = feval(func,'f',sstar,xstar,estar,params{:});
    [g,gx,gxx,gs] = feval(func,'g',sstar,xstar,estar,params{:});
  end
end

% RESHAPE TO ENSURE COMFORMABILIITY
fs  = reshape(fs , 1,ds);
fx  = reshape(fx , 1,dx);
fss = reshape(fss,ds,ds);
fxs = reshape(fxs,dx,ds);
fxx = reshape(fxx,dx,dx);
g   = reshape(g  ,ds, 1);
gx  = reshape(gx ,ds,dx);
gs  = reshape(gs ,ds,ds);
gxx = reshape(delta*pstar*reshape(gxx,ds,dx*dx),dx,dx);
if foc
  gxs = reshape(delta*pstar*reshape(gxs,ds,dx*ds),dx,ds);
  gss = reshape(delta*pstar*reshape(gss,ds,ds*ds),ds,ds);
else
  gxs = 0;
  gss = 0;
end

if qz % Use QZ decomposition
  if foc
    A=[ eye(ds)          zeros(ds,dx+ds); 
        zeros(dx,ds+dx) -delta*gx'      ;
        zeros(ds,ds+dx)  delta*gs'      ];
    B=[ gs       gx        zeros(ds,ds);
        fxs+gxs  fxx+gxx   zeros(dx,ds);
       -fss-gss -fxs'-gxs' eye(ds)     ];
  else 
    A=[ eye(ds)          zeros(ds,dx+ds); 
        zeros(dx,ds+dx) -delta*gx'      ;
        zeros(ds,ds+dx)  delta*gs'      ];
    B=[ gs   gx   zeros(ds,ds);
        fxs  fxx  zeros(dx,ds);
       -fss -fxs' eye(ds)     ];
  end
  [S,T,Q,Z]=qzordered(A,B);
  C=real(Z(ds+1:end,1:ds)/Z(1:ds,1:ds));
  Gamma=C(1:dx,:);
  Delta=C(dx+1:end,:);
  G=gs+gx*Gamma;
else % Use Riccati equation function iteration
  tol=1.e-8;
  % Starting values
  if cond(gx'*gx)>1.e10
    Delta = eye(ds);
  else
    Delta = -((delta*gx'*gx)...
             \(((delta*gs'*gx)*(delta*gx'*gs))+delta*gs'*gs-fxx'));
  end
  % Iterate on Riccati equation
  for it=1:200
    Deltaold = Delta;
    temp  = delta*gx'*Delta*gs+fxs;
    Gamma = -((delta*gx'*Delta*gx+fxx)\temp);
    Delta = fss+delta*gs'*Delta*gs+temp'*Gamma;
    if norm(Delta-Deltaold,inf)<tol, break, end;
  end
  G=gs+gx*Gamma;
end

n = size(s,1);
s = s-sstar(ones(1,n),:);

% COMPUTE CONTROL, SHADOW PRICE AND VALUE FUNCTIONS AT REQUESTED EVALUATION POINTS
xlq = xstar(ones(1,n),:) + s*Gamma';
plq = pstar(ones(1,n),:) + s*Delta';
vlq = feval(func,'f',sstar,xstar,estar,params{:})/(1-delta) ...
            + s*pstar' + 0.5*sum(s.*(s*Delta'),2); 

