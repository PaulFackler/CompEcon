% DEMOPT04 Demonstrates Quasi-Newton maximization
% Maximizing banana function f(x,y)=-100*(y-x*x)^2-(1-x)^2, starting at [0;1]. 
function demopt04
close all

disp(' ')
disp('DEMOPT04 Demonstrates Quasi-Newton maximization')
warning off

n = [20 20];
xmin = [-0.2 -0.2];
xmax = [ 1.2  1.2];
[x,xcoord] = nodeunif(n,xmin,xmax);
[x1,x2] = meshgrid(xcoord{1},xcoord{2}) ;

y = banana(x');
y = reshape(y,n(1),n(2))';
conts = -exp(0.25:0.5:20);

optset('qnewton','maxit',1);
optset('qnewton','ShowIters',1);

'Steepest Ascent Maximization'
optset('qnewton','SearchMeth',1);
k = 100;
x0 = [0;1];
x=x0;
A=[];
for i=1:k
  xx1(:,i) = x;
  [x,A] = qnewton('banana',x,A);
  if norm(x-[1;1])>sqrt(eps), iters1=i; end
end

'DFP Maximization'
optset('qnewton','SearchMeth',2);
x=x0;
A=[];
for i=1:k
  xx2(:,i) = x;
  [x,A] = qnewton('banana',x,A);
  if norm(x-[1;1])>sqrt(eps), iters2=i; end
end

'BFGS Maximization'
optset('qnewton','SearchMeth',3);
x=x0;
A=[];
for i=1:k
  xx3(:,i) = x;
  [x,A] = qnewton('banana',x,A);
  if norm(x-[1;1])>sqrt(eps), iters3=i; end
end

figure(1)
plot(1,1,'o')
hold on
plot(xx1(1,:),xx1(2,:))
plot(xx1(1,:),xx1(2,:),'*')
contour(xcoord{1},xcoord{2},y,conts,'k:')
hold off
axis square
title('Steepest Ascent Maximization of Banana Function')
h=xlabel('x_1');
set(h,'VerticalAlignment','cap')
h=ylabel('x_2');
set(h,'VerticalAlignment','bottom')
xlim([xmin(1) xmax(1)])
ylim([xmin(2) xmax(2)])
set(gca,'ytick',[0 0.5 1])

figure(2)
plot(1,1,'o')
hold on
plot(xx2(1,:),xx2(2,:))
plot(xx2(1,:),xx2(2,:),'*')
contour(xcoord{1},xcoord{2},y,conts,'k:')
hold off
axis square
title('DFP Quasi-Newton Maximization of Banana Function')
h=xlabel('x_1');
set(h,'VerticalAlignment','cap')
h=ylabel('x_2');
set(h,'VerticalAlignment','bottom')
xlim([xmin(1) xmax(1)])
ylim([xmin(2) xmax(2)])
set(gca,'ytick',[0 0.5 1])

figure(3)
plot(1,1,'o')
hold on
plot(xx3(1,:),xx3(2,:))
plot(xx3(1,:),xx3(2,:),'*')
contour(xcoord{1},xcoord{2},y,conts,'k:')
hold off
axis square
title('BFGS Quasi-Newton Maximization of Banana Function')
h=xlabel('x_1');
set(h,'VerticalAlignment','cap')
h=ylabel('x_2');
set(h,'VerticalAlignment','bottom')
xlim([xmin(1) xmax(1)])
ylim([xmin(2) xmax(2)])
set(gca,'ytick',[0 0.5 1])
s={'Iterations';
   ['Steepest: ' num2str(iters1+1)];...
   ['DFP:      ' num2str(iters2+1)];...
   ['BFGS:     ' num2str(iters3+1)]};
disp(s)

prtfigs(mfilename,'Steepest Ascent Maximization of Banana Function',1)
prtfigs(mfilename,'BFGS Quasi-Newton Maximization of Banana Function',3)
