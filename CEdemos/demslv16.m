% DEMSLV16 Billups Hard NCP
disp('DEMSLV16 Billups Hard NCP')
close all

n=10; 
xx=zeros(n,1);
optset('newton','maxit',1);
optset('newton','maxsteps',50);
warning off;
for i=2:n
   xx(i)=newton('billupss',xx(i-1)); 
end
warning on;
optset('newton','maxit',100);

close all
figure(1);
subplot(1,2,1);
x=nodeunif(500,-0.5,2.5);
plot(x,billupss(x),'k-',x,billupsm(x),'k--',x,zeros(size(x)),'k:');
xlabel('x')   
axis([-0.5 2.5 -1.0 1.5])
axis square
title('A Difficult NCP')

subplot(1,2,2);
x=nodeunif(500,-0.03,0.03);
plot(x,billupss(x),'k-',x,billupsm(x),'k--',x,zeros(size(x)),'k:');
xlabel('x') 
title('A Difficult NCP Magnified')
axis([-.03 .03 -.01 .05])
axis square

optset('newton','defaults')

prtfigs(mfilename,'Billups'' Hard NCP',1)
