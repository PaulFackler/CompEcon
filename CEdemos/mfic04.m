% MFIC04 Model function file for capacity choice example
% See DEMIC04 for usage
function out=mfic04(flag,s,P,C,delta,rho)
switch flag
case 'f'
  out=P*log(s+1);
case 'mu'
  out=-delta*s;
case 'sigma'
  out=[];
case 'rho'
  out=rho+zeros(size(s,1),1);
case 'R+'
  out=[C*(s(1)-s(2)) C -C 0];
case 'R-'
  out=zeros(1,4);
end
