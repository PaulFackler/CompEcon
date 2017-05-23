% REMSIMUL  Monte Carlo simulation of discrete time controlled Markov process
% USAGE
%   st = remsimul(model,s0,nper,sres,xres,)
% INPUTS
%   model   : name of model structure
%   fspace  : name of projection space structure
%   s0      : k by 1 vector of initial states
%   nper    : number of simulated time periods
%   sres    : coordinates of the evaluation grid (from remsolve)
%   xres    : optimal control function values 
%               at grid defined by sres (from remsolve)
% OUTPUT
%   st      : k by nyrs vector of simulated states
%
% USES: DISCRAND

% Copyright (c) 1997-2000, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function [ssim,xsim] = remsimul(model,s0,nper,sres,xres)

% DUMMY SHOCKS/WEIGHTS IF MODEL DETERMINISTIC  
  if isfield(model,'e'), e=model.e; else, e=0; end;
  if isfield(model,'w'), w=model.w; else, w=1; end;
  
  nrep = size(s0,1);
  ds   = size(s0,2);
  if ds>1
    st   = gridmake(sres);
    dx   = ds*length(xres(:))/length(st(:));
  else
    dx   = length(xres(:))/length(sres(:));
  end
  ssim = zeros(nrep,ds,nper+1);  
  xsim = zeros(nrep,dx,nper+1);
  
  nx = prod(size(xres))/dx;
  xres = reshape(xres,nx,dx);
  for t=1:nper+1
    xx = minterp(sres,xres,s0);
    ssim(:,:,t) = s0;      
    xsim(:,:,t) = xx;
    ee = e(discrand(nrep,w),:);
    s0 = feval(model.func,'g',s0,xx,[],ee,model.params{:});
  end
 
  ssim = ssim(:,:,1:nper+1); 
  xsim = xsim(:,:,1:nper+1);
  if ds==1; ssim=squeeze(ssim(:,1,:)); end 
  if dx==1; xsim=squeeze(xsim(:,1,:)); end