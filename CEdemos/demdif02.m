% DEMDIF02 Demonstrates errors in finite difference differentiation of e(x)
function demdif02
close all

disp(' ')
disp('DEMDIF02 Demonstrates errors in finite difference differentiation of e(x)')

n=50;
df=zeros(n,1);
h=ones(n,1);
x=1;
for i=1:n;
  if i>1;h(i)=h(i-1)/2;end;
  df(i)=(exp(x+h(i))-exp(x))./h(i);
end;
figure(1)
plot(log10(h),log10(abs(df-exp(x))),'k',log10(sqrt(eps)),-10,'k*');
title('Errors in 1-Sided Numerical Derivatives')
hh=xlabel('log_{10}(h)');
pp=get(hh,'position'); pp(2)=pp(2)+0.25;
set(hh,'position',pp);
ylabel('log_{10}(e)') 

for i=1:n;
  if i>1;h(i)=h(i-1)/2;end;
  df(i)=(exp(x+h(i))-exp(x-h(i)))./(2*h(i));
end;
figure(2)
plot(log10(h),log10(abs(df-exp(x))),'k',log10(eps.^(1/3)),-12,'k*');
title('Errors in 2-Sided Numerical Derivatives')
hh=xlabel('log_{10}(h)');
pp=get(hh,'position'); pp(2)=pp(2)+0.25;
set(hh,'position',pp);
ylabel('log_{10}(e)') 

prtfigs(mfilename,'Errors in Finite Difference Numerical Derivatives',[1 2])
