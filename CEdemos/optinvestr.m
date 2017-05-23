% OPTINVESTR Residual function for optimal investment problem
function [e,c]=optinvestr(Sstar,I,B,B1,phi11,n) 

  c=(B-Sstar*B1)\[zeros(n-1,1);Sstar-I];
  e=Sstar-phi11*c;
