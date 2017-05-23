function [fx,J]=billupss(x)
% BILLUPSM Semismooth reformulation of Billups' function
[fx,J] = billups(x);
if nargout<2
   fx = smooth(x,0,inf,fx,J);
else
   [fx,J] = smooth(x,0,inf,fx,J);
end

