% SCSTEADYSTATE Computes deterministic steady state values
% USAGE
%  [s,lambda,x]=scsteadystate(model,s0,lambda0);
% INPUTS
%  model : a model structure as defined in SCSOLVE
%  s0    : an initial guess for the steady state value of the state variable
%  lambda0 : an initial guess for the steady state value of Vs
% OUTPUTS
%  s      : steady state value of the state variable
%  lambda : steady state value of Vs (co-state variable)
%  x      : steady state value of the control
%
% See SCSOLVE for the definition of the model structure
function [s,lambda,x]=scsteadystate(model,s0,lambda0)

func = model.func;     % model functions
params=model.params;   % other parameters

d=length(s0);
z=broyden(@scroots,[s0(:);lambda0(:)],func,d,params);
s=z(1:d);
lambda=z(d+1:end);
x=feval(func,'x',s',[],lambda',params{:})';

function res=scroots(z,func,d,params)
s=z(1:d);
lambda=z(d+1:end);
x=feval(func,'x',s',[],lambda',params{:});
g=feval(func,'g',s',x,lambda',params{:});
df=fjac(func,[2],'f',s',x,lambda',params{:});
dg=fjac(func,[2],'g',s',x,lambda',params{:});
rho=feval(func,'rho',s',x,lambda',params{:});
res=[g';(rho*lambda'-lambda'*dg-df)'];
return
