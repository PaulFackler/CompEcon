% REMSTST  Computes invariant distribution for rational expectations models
% USAGE
%   [smid,pi,xx] = remstst(model,nsmooth,nbin,sres,x);
% INPUTS
%   model   : name of model structure
%   nsmooth : positive integer, smoothing parameter
%   nbin    : number of bins in distribution histogram
%   sres    : coordinates of the evaluation grid (from remsolve)
%   x       : responses at states sres (from remsolve)
% OUTPUTS
%   smid    : midpoints of distribution histogram
%   pi      : associated invariant probabilities
%   xx      : responses at states smid
%
% Setable options (use OPTSET):
%   tol     : convergence tolerance
%   maxit   : maximum number of iterations

% Copyright (c) 1997-2002, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function [smid,pi,xx] = remstst(model,nsmooth,nbin,s,x)

maxit = optget('remstst','maxit',100);
tol   = optget('remstst','tol',sqrt(eps));
if ~isfield(model,'e'); model.e=0; end;
if ~isfield(model,'w'); model.w=1; end;

smin = min(s);
smax = max(s);

m = length(model.w);
swid = (smax-smin)/nbin;
smid = nodeunif(nbin,smin+swid/2,smax-swid/2); 
ss = smid; qq = ones(size(smid));
for t=1:nsmooth
  ns = length(ss);
  xx = minterp(s,x,ss);
  ss = ss(:,ones(1,m));
  xx = xx(:,ones(1,m));
  qq = qq(:,ones(1,m));
  ee = model.e(:,ones(1,ns))';
  ww = model.w(:,ones(1,ns))';
  qq = qq.*ww;
  ss = ss(:);
  xx = xx(:);
  ee = ee(:);
  qq = qq(:);
  ss = feval(model.func,'g',ss,xx,[],ee,model.params{:});
end
nobs = m^nsmooth;
ss = reshape(ss,nbin,nobs);
qq = reshape(qq,nbin,nobs);
P = zeros(nbin,nbin);
for i=1:nbin
  j = ceil((ss(i,:)-smin)/swid);
  j = min(j,nbin); j = max(j,1);
  for k=1:nobs;
    P(i,j(k)) = P(i,j(k)) +  qq(i,k);
  end
end

xx = minterp(s,x,smid);
pi = markov(P);
