% Model function file for commodity storage example
function out=pbvp02(flag,t,x,dx,r,C,eta,S0);
  switch flag
  case 'r'
    out=dx-[r*x(:,1)+C  -x(:,1).^(-eta)];
  case 'b'
    out=x(:,2)-[S0;0];
  end
  
