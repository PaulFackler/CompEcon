% BANANA The "banana" or Rosencrantz function
function [y,dy]=banana(x)
 y  = -100*(x(2,:)-x(1,:).^2).^2-(1-x(1,:)).^2;
 if nargout>1
 dy = [2*(1-x(1,:)) + 400*(x(2,:)-x(1,:).^2).*x(1,:); ...
       -200*(x(2,:)-x(1,:).^2)];
 end
