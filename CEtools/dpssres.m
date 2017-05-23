% DPSSRES Residual function for computing steady state in DP Models
% See DPSS

% Copyright (c) 1997-2002,  Paul L. Fackler & Mario J. Miranda 
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function r=ressres(theta,func,params,e,delta,ds,dx);

   se=theta(1:ds)';
   xe=theta(ds+1:ds+dx)';
   lambdae=theta(ds+dx+1:end)';
   
  if nargout(func)<4    % only x derivatives available
    [f,fx] = feval(func,'f',se,xe,e,params{:});
    [g,gx] = feval(func,'g',se,xe,e,params{:});
    fs  = fjac(func,2,'f',se,xe,e,params{:});
    gs  = fjac(func,2,'g',se,xe,e,params{:});
  else                  % s and x derivatives available
    [f,fx,fxx,fs] = feval(func,'f',se,xe,e,params{:});
    [g,gx,gxx,gs] = feval(func,'g',se,xe,e,params{:});
  end
  
  fs  = reshape(fs, 1,ds);
  fx  = reshape(fx, 1,dx);
  gx  = reshape(gx,ds,dx);
  gs  = reshape(gs,ds,ds);
  
  
  rx=fx+delta*lambdae*gx;
  [a,b]=feval(func,'b',se,xe,e,params{:});
  rx=min(max(rx,a-xe),b-xe);
  r=[se-g rx lambdae-(fs+delta*lambdae*gs)]';
