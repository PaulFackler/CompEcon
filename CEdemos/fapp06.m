% FAPP06 A 1-D function to be approximated
function [y,d,s] = fapp06(x);
y = exp(-x);
d = -exp(-x);
s = exp(-x);
