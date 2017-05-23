% MFRE03 Model function file for government price control
function out = mfre02(flag,s,x,z,e,snext,xnext,delta,xmax,pstar,gamma,c0,c1);
switch flag
case 'f'; 
   out(:,1) = pstar - (s-x(:,1)).^(-gamma);                % f
   out(:,2) = delta*z - (c0+c1*x(:,2));  
case 'x'
  out = max([s-pstar^(-1/gamma) (delta*z-c0)/c1],0);
case 'g'; 
  out =  x(:,1) + x(:,2).*e;                                % g
case 'h'
  out=(snext-xnext(:,1)).^(-gamma).*e;
case 'b'; % BOUND FUNCTION
   out1 = zeros(n,2);                                % xl
   out2 = [xmax*ones(n,1) inf*ones(n,1)];            % xu
end 