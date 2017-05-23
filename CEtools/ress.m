% RESS Solves for the steady state in rational expectations models
% USAGE
%   [se,xe,ze]=ress(model,s0,x0,z0)
% INPUTS
%   model    : structure variable (see RESOLVE)
%   s0,x0,z0 : initial values of variables
% OUTPUTS
%   se,xe,ze : equilibrium values of variables

% Copyright (c) 1997-2002,  Paul L. Fackler & Mario J. Miranda 
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function [se,xe,ze]=ress(model,s0,x0,z0)

 func=model.func;
 params=model.params;
 e=model.w'*model.e;
 n=length(s0);
 m=length(x0);
 if nargin>3, p=length(z0); else, p=[]; end

 optset('broyden','initi',0);
 optset('broyden','maxit',100);

 if model.explicit & 0
   x0=feval(func,'x',s0,[],z0,[],[],[],params{:});
   theta=[s0(:);z0(:)];
   theta=real(broyden('ressres',theta,func,params,e,n,m,p,1));
   se=theta(1:n)';
   ze=theta(n+1:end)';
   xe=feval(func,'x',se,[],ze,[],[],[],params{:});
 elseif 0
   theta=[s0(:);x0(:);z0(:)];
   theta=real(broyden('ressres',theta,func,params,e,n,m,p,0));
   se=theta(1:n)';
   xe=theta(n+1:n+m)';
   ze=theta(n+m+1:end)';
 else
   theta=[s0(:);x0(:)];
   theta=real(broyden('ressres',theta,func,params,e,n,m,p,0));
   se=theta(1:n)';
   xe=theta(n+1:n+m)';
   ze=feval(func,'h',se,xe,[],e,se,xe,params{:});
 end
 
 optset('broyden','defaults');
