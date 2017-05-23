% FAPP10 Residual function for function inverse collocation example
function resid=fapp10(c,x,Phi)
resid = exp(Phi*c)-x;
