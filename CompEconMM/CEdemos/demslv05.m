%% DEMSLV05 Newton and Broyden Paths
%
% Plots path taken by Newton and Broyden iterates, starting from initial
% value (1.0.,0.5), in solving bivariate rootfining problem
% x2*exp(x1)-2*x2=0 and x1*x2-x2^3=0.

function demslv05

% Preliminary tasks
demosetup(mfilename)
warning off

% Set Initial Value and Max Iterations
xinit = [1.0;0.5];
niter = 10;

% Compute Newton Path
optset('newton','showiters',0)
optset('newton','maxit',1)
xn = zeros(2,niter);
xn(:,1) = xinit;
for i=2:niter
  xn(:,i) = newton(@f,xn(:,i-1));
end

% Compute Broyden Path
optset('broyden','showiters',0)
optset('broyden','maxit',1)
xb = zeros(2,niter);
xb(:,1) = xinit;
fjacinv = inv(fdjac(@f,xinit));
for i=2:60
  optset('broyden','initb',fjacinv)
  [xb(:,i),~,fjacinv] = broyden(@f,xb(:,i-1));
end

% Set figure parameters
figure
x1min = 0.3;
x1max = 1.1;
x2min = 0.4;
x2max = 1.2;
xticks = [];
yticks = [];

% Generate data for contour plot
n = 100;
x1 = nodeunif(n,x1min,x1max);
x2 = nodeunif(n,x2min,x2max);
z1 = zeros(n,n);
z2 = zeros(n,n);
for i=1:n
  for j=1:n
    z = f([x1(j);x2(i)]);
    z1(i,j) = z(1);
    z2(i,j) = z(2);
  end
end

% Subplot for Newton Path
subplot(1,2,1)
hold on
contour(x1,x2,z1,[0 0],'k','Linewidth',2)
contour(x1,x2,z2,[0 0],'k','Linewidth',2)
for i=1:niter-1
  plot([xn(1,i) xn(1,i+1)],[xn(2,i) xn(2,i+1)],'b-')
  plotbullet(xn(1,i),xn(2,i),18,'r')
end
plotbullet(xn(1,niter),xn(2,niter),18,'r')
plottext(0.70,0.45,'$f_1=0$',[],[],8)
plottext(0.47,0.65,'$f_2=0$',[],[],8)
set(gca,'xtick',xticks)
set(gca,'ytick',yticks)
axis('square')
title('Newton''s Method')
xlabel('$x_1$')
ylabel('$x_2$')

% Subplot for Broyden Path
subplot(1,2,2)
hold on
contour(x1,x2,z1,[0 0],'k','Linewidth',2)
contour(x1,x2,z2,[0 0],'k','Linewidth',2)
for i=1:niter-1
  plot([xb(1,i) xb(1,i+1)],[xb(2,i) xb(2,i+1)],'b-')
  plotbullet(xb(1,i),xb(2,i),18,'r')
end
plotbullet(xb(1,niter),xb(2,niter),18,'r')
plottext(0.70,0.45,'$f_1=0$',[],[],8)
plottext(0.47,0.65,'$f_2=0$',[],[],8)
set(gca,'xtick',xticks)
set(gca,'ytick',yticks)
axis('square')
title('Broyden''s Method')
xlabel('$x_1$')
ylabel('$x_2$')

% SAVE FIGURES
printfigures(mfilename,1)


function [fval,fjac] = f(x)
fval = [x(2)*exp(x(1))-2*x(2);x(1)*x(2)-x(2)^3];
fjac = [x(2)*exp(x(1)) exp(x(1))-2; x(2) x(1)-3*x(2)^2];