% LCPSTEP - Newton step for Array Linear Complementarity Problem

% Copyright (c) 1997-2002, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function [F,dx] = lcpstep(method,x,xl,xu,F,Fx)
if method(1)=='m'
  if nargout==1
    F = min(max(F,xl-x),xu-x);
  else  
    F = min(max(F,xl-x),xu-x);
    dx = -arrayinvb(F,Fx,x,xl,xu);
    dx = min(max(dx,xl-x),xu-x);
  end
else
  if nargout==1
    F = arrayss(x,xl,xu,F);
  else  
    [F,Fx] = arrayss(x,xl,xu,F,Fx);
    dx = -arrayinv(F,Fx);
    dx = min(max(dx,xl-x),xu-x);
  end
end


