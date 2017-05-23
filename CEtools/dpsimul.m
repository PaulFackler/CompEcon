% DPSIMUL  Monte Carlo simulation of discrete time controlled Markov process
% USAGE
%   [ssim,xsim] = dpsimul(model,ss,nper,sres,xres)
% INPUTS
%   model   : model structure variable
%   ss      : k by d vector of initial states
%   nper    : number of simulated time periods
%   sres    : coordinates of the evaluation grid (from dpsolve)
%   xres    : optimal control function values 
%               at grid defined by sres
% OUTPUTS
%   ssim    : k by d by nper+1 vector of simulated states       
%   xsim    : k by d by nper+1 vector of simulated actions
% For finite horizon problems, xsim is k by d by nper
% If d=1, ssim and xsim are k by nper+1
%
% USES: DISCRAND

% Copyright (c) 1997-2002, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function [ssim,xsim] = dpsimul(model,ss,nper,sres,xres)

% DUMMY SHOCKS/WEIGHTS IF MODEL DETERMINISTIC  
  if isfield(model,'horizon'), T=model.horizon; else, T=inf; end;
  if isfield(model,'e'), e=model.e; else, e=0; end;
  if isfield(model,'w'), w=model.w; else, w=1; end;

  func=model.func;
  params=model.params;
  
  nper = min(nper,T);
  nrep = size(ss,1);
  ds   = size(ss,2);
  st   = gridmake(sres);
  if T==inf
    dx = ds*length(xres(:))/length(st(:));
  else
    dx = ds*length(xres(:))/(T*length(st(:)));
  end
  ssim = zeros(nrep,ds,nper+1);  
  xsim = zeros(nrep,dx,nper+1);  

  if T<inf
    nx = numel(xres)/(dx*T);
    xres = reshape(xres,nx,T,dx);
    ssim(:,:,1) = ss; 
    for t=1:nper
      xx = minterp(sres,xres(:,t,:),ss); 
      xsim(:,:,t) = xx; 
      ee = e(discrand(nrep,w),:);
      ss = feval(func,'g',ss,xx,ee,params{:});
      ssim(:,:,t+1) = ss; 
    end
  else
    nx = numel(xres)/dx;
    xres = reshape(xres,nx,dx);
    for t=1:nper+1
      xx = minterp(sres,xres,ss);
      ssim(:,:,t) = ss;      
      xsim(:,:,t) = xx;
      ee = e(discrand(nrep,w),:);
      ss = feval(func,'g',ss,xx,ee,params{:});
    end   
 end
 
 if T<inf
   xsim(:,:,end) = [];
 end

 if ds==1; ssim=squeeze(ssim); end 
 if dx==1; xsim=squeeze(xsim); end