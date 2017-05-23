% Model function file BVP example
function out=pbvp01(flag,t,x,dx,A);
  switch flag
  case 'r'
    out=dx-x*A;
  case 'b'
    out=[x(1,1)-1;x(2,2)-1];
  end
  
