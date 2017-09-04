%% DEMMATH01 Taylor Approximations
%
% Illustrate uni- and bi-variate Taylor approximations.

% Preliminary tasks
demosetup(mfilename)


%% Univariate Taylor Approximation

% Compute first- and second-order approximations
x = nodeunif(100,-1,1);
y = (x+1).*exp(2*x);
y1 = 1+3*x;
y2 = 1+3*x+4*x.^2;

% Plot first- and second-order approximations
figure
hold on
plot(x,y ,'k')
plot(x,y1,'b')
plot(x,y2,'r')
legend('Function','1st Order Approximation','2nd Order Approximation')
title('Taylor Approximations for Univariate Fuction')
xlabel('$x$');
ylabel('$y$');

%% Bivariate Taylor Approximation

% Compute first- and second-order approximations
[x1,x2] = ndgrid(0:0.1:2,-1:0.01:1);
f  = (x1.^2).*exp(-x2);
f1 = 2*x1 - x2 - 1;
f2 = x1.^2 - 2*x1.*x2 + 0.5*x2.^2 + x2;

% Plot first-order approximation error
figure
mesh(x1,x2,f1-f) 
title('First-Order Approximation Error')
xlabel('$x_1$')
ylabel('$x_2$')
zlabel('$y$') 

% Plot second-order approximation error
figure
mesh(x1,x2,f2-f)
title('Second-Order Approximation Error')
xlabel('$x_1$')
ylabel('$x_2$')
zlabel('$y$')

%% SAVE FIGURES
printfigures(mfilename)