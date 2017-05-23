% MFDP13 Function file for production-inventory model
function [out1,out2,out3] = mfdp13(flag,s,x,e,c,k,pbar,rho);
n = size(s,1);
ds = 2;
dx = 2;
switch flag
case 'b'
  out1 = zeros(size(s));
  out2 = inf*ones(size(s));
case 'f'
  out1 = zeros(n,1); 
  out2 = zeros(n,dx);
  out3 = zeros(n,dx,dx);
  out1 = s(:,2).*(s(:,1)+x(:,1)-x(:,2)) ...
        - (c(1)+0.5*c(2)*x(:,1)).*x(:,1) ...
        - (k(1)+0.5*k(2)*x(:,2)).*x(:,2);
  out2(:,1) =  s(:,2) - (c(1)+c(2)*x(:,1));
  out2(:,2) = -s(:,2) - (k(1)+k(2)*x(:,2));
  out3(:,1,1) = -c(2)*ones(n,1);
  out3(:,2,2) = -k(2)*ones(n,1);
case 'g'
  out1 = zeros(n,ds);
  out2 = zeros(n,ds,dx);
  out3 = zeros(n,ds,dx,dx);
  out1 = [x(:,2) pbar+rho*(s(:,2)-pbar)+e];
  out2(:,1,2) = ones(n,1);
end