% FDIF03 ODE file for commercial fishery example
% Defines the differential equation used in demdif03.m

function dx=fdif03(t,x,flag,beta,f,del);  
  s=x(1,:);
  k=x(2,:);
  temp=1+beta*s.*k;
  ds=(1-s).*s-s.*k./temp;
  dk=del*(s./(temp.^2)-2*f);
  dx=[ds;dk];
