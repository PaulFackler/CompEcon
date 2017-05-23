% DPSS Solves for the steady state in dynamic programming models
% USAGE
%   [se,xe,lambdae]=dpss(model,s0,x0,lambda0)
% INPUTS
%   model : structure variable (see DPSOLVE)
%   s0,x0,lambda0 : initial values of variables 
%                   (state, control and shadow price)
% OUTPUTS
%   sstar,xstar,lambdastar : steady state values of variables

% Copyright (c) 1997-2002,  Paul L. Fackler & Mario J. Miranda 
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function [sstar,xstar,lambdastar]=dpss(model,s0,x0,lambda0)
  if nargin<5
    if isfield(model,'w'), e=model.w'*model.e; 
    else e=[];
    end
  end

 func=model.func;
 params=model.params;

 delta=model.discount;
 n=length(s0);
 m=length(x0);

 optset('broyden','initi',0);
 optset('broyden','maxit',100);

 theta=[s0(:);x0(:);lambda0(:)];
 theta=real(broyden('dpssres',theta,func,params,e,delta,n,m));
 sstar=theta(1:n)';
 xstar=theta(n+1:n+m)';
 lambdastar=theta(n+m+1:end)';