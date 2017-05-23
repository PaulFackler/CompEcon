% FAPP05 A set of functions to demonstrate approximation methods
function y = fapp05(x,i);
switch i
   case 1, y = 1 + x + 2*x.^2 - 3*x.^3;
   case 2, y = exp(-x);
   case 3, y = 1./(1+25*x.^2);
   case 4, y = sqrt(abs(x));
end


