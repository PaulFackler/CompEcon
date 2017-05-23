% SCRES Computes residual functions for stochastic control problems
% USAGE
%   r=scres(model,fspace,s,cv,cx)
%
% Univariate only. See SCSolve for way to make multivariate

% Copyright (c) 1997-2000, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function r=scres(model,fspace,s,cv,cx)

  scfile=model.func;

  V=funeval(cv,fspace,s,[0;1;2]);
  r=feval(scfile,'rho',s,[],[],model.params{:}).*V(:,:,1);
  r=r-0.5*feval(scfile,'sigma',s,[],[],model.params{:}).^2.*V(:,:,3);
  if nargin<6 | isempty(cx)
    x=feval(scfile,'x',s,[],V(:,:,2),model.params{:});
  else
    x=funeval(cx,fspace,s); 
  end
  f=feval(scfile,'f',s,x,[],model.params{:});
  g=feval(scfile,'g',s,x,[],model.params{:});
  r=r-g.*V(:,:,2)-f; 