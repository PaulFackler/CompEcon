% DPSTST  Computes invariant distribution for continuous-state/action controlled dynamic program
% USAGE
%   [smid,p,xx] = dpstst(model,nsmooth,nnodes,scoord,x);
% INPUTS
%   model   : name of model structure
%   nsmooth : positive integer, smoothing parameter
%   nnodes  : number of nodal values in distribution approximation
%   scoord  : coordinates of the evaluation grid (from dpsolve)
%   x       : optimal control function values 
%               at grid defined by scoord (from dpsolve)
% OUTPUTS
%   smid    : midpoints of distribution histogram
%   p       : associated invariant probabilities
%   xx      : control values at smid
%
% Estimates of the steady state density of the state variable
% can be obtained as 
%    p/binwidth, where binwidth=(scoord(end)-scoord(1))/nnodes
%
% NOTE: only implemented for 1-dimensional state space

% Copyright (c) 1997-2000, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function [smid,p,xx] = dpstst(model,nsmooth,nnodes,scoord,x)

d=size(scoord,2);
if d>1, error('Implemented only for 1-D states'); end

e=model.e;
w=model.w;
func=model.func;
params=model.params;

smin=scoord(1);
smax=scoord(end);
swid = (smax-smin)./nnodes;
smid=linspace(smin+swid/2,smax-swid/2,nnodes)';

m = length(w);
ss=gridmake(smid);
q = ones(size(ss));
for t=1:nsmooth
  xx = minterp(scoord,x,ss);
  [ss,xx,ee]=gridmake([ss xx],e);
  ss = feval(func,'g',ss,xx,ee,params{:});
  q = q*w'; q=q(:);
end
nobs = m^nsmooth;
ss = reshape(ss,nnodes,nobs);
q = reshape(q,nnodes,nobs);

p = zeros(nnodes,nnodes);
for j=1:nobs
  i = ceil((ss(:,j)-smin)/swid);
  i = min(max(i,1),nnodes);
  ind=(1:nnodes)'+(i-1)*nnodes;
  p(ind) = p(ind) +  q(i,j);
end
p = markov(p);

if nargout>2, xx = minterp(scoord,x,smid); end
