% DEMAPP00 Creates Table 6.1

runge = inline(' 1./(1+25*x.^2)','x');

n = 10:10:50;
a = -5;
b =  5;
x = linspace(-5,5,501)';

for i=1:length(n)
  fspace=fundefn('cheb',n(i),a,b);
  nodes=linspace(a,b,n(i))';
  c=funbas(fspace,nodes)\runge(nodes);
  yhat=funeval(c,fspace,x);
  fprintf('%2i   %5.2f\n',n(i),log10(norm(runge(x)-yhat,inf)))
end

