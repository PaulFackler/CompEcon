function [fx,J]=billupsm(x)
% BILLUPSM Minmax reformulation of Billups' function
[fx,J] = billups(x);
if nargout<2
   fx = minmax(x,0,inf,fx,J);
else
   [fx,J] = minmax(x,0,inf,fx,J);
end

