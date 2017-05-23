% DEMOPT03 Demonstrates Nelder-Mead simplex method
% Creates and plays a movie of Nelder-Meade simplex iterations when
% maximizing banana function f(x,y)=-100*(y-x*x)^2-(1-x)^2, starting at [0;1]. 
% To view the movie again use
%   M=demopt03; movie(M);
function M=demopt03
close all

disp(' ')
disp('DEMOPT03 Demonstrates Nelder-Mead simplex method')

n = [20 20];
xmin = [-0.2 -0.2];
xmax = [ 1.2  1.2];
[x,xcoord] = nodeunif(n,xmin,xmax);
[x1,x2] = meshgrid(xcoord{1},xcoord{2}) ;

y = banana(x');
y = reshape(y,n(1),n(2))';
conts = -exp(0.25:0.5:20);

figure(1)
contour(xcoord{1},xcoord{2},y,conts,'k:')
xlabel('x_1'),ylabel('x_2')
title('Nelder-Mead Maximizes the Banana Function')

optset('neldmead','maxit',1);
k = 50;
x = [1;0];
warning off
[xx,S] = neldmead('banana',x);
warning on
hold on
hp = patch(S(1,:),S(2,:),[0.5 0.5 0.5]);
M  = moviein(k);
for i=1:k
  xvec(:,i) = x;
  warning off
  [x,S] = neldmead('banana',x,S);
  warning on
  set(hp,'xdata',S(1,:)','ydata',S(2,:)');
  M(:,i) = getframe;
end
hold off
optset('neldmead','defaults');

%for i=1:size(M,2),movie(M(:,i),0);end

figure(2)
plot(1,1,'o')
hold on
plot(xvec(1,:),xvec(2,:))
plot(xvec(1,:),xvec(2,:),'*')
contour(xcoord{1},xcoord{2},y,conts,'k:')
hold off
axis square
title('Nelder-Mead Maximization of Banana Function')
h=xlabel('x_1');
set(h,'VerticalAlignment','cap')
h=ylabel('x_2');
set(h,'VerticalAlignment','bottom')
axis([-.2 1.2 -.2 1.2])
set(gca,'ytick',[0 0.5 1])

prtfigs(mfilename,'Nelder-Mead Maximization of Banana Function',2)
