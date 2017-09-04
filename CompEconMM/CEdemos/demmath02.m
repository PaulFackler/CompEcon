%% DEMMATH02 Function Inner Products, Norms & Metrics
%
% Illustrates computation of inner product, norm, metric for a cunction. 

% Preliminary tasks
demosetup(mfilename)

% Compute Inner Product and Angle
l=-1; u=1;
f = @(x) 2*x.^2-1;
g = @(x) 4*x.^3-3*x;
fg = integral(@(x)f(x).*g(x),l,u);
ff = integral(@(x)f(x).*f(x),l,u);
gg = integral(@(x)g(x).*g(x),l,u);
angle = acosd(fg/sqrt(ff*gg));

% Compute Function Norm
l=0; u=2;
f = @(x) x.^2-1;
p = 1; integral(@(x)abs(f(x)).^p,l,u)^(1/p);
p = 2; integral(@(x)abs(f(x)).^p,l,u)^(1/p);

% Compute Function Metrics
l=0; u=1;
f = @(x) 5+5*x.^2;
g = @(x) 4+10*x-5*x.^2;
p = 1; integral(@(x)abs(f(x)-g(x)).^p,l,u)^(1/p);
p = 2; integral(@(x)abs(f(x)-g(x)).^p,l,u)^(1/p);

% Illustrate Function metrics
x = nodeunif(200,l,u);
figure
hold on
xlabel('$x$')
ylabel('$|f(x)-g(x)|$')
set(gca,'XTick',[0 1 2])
set(gca,'YTick',[0 1 2 3])
plot(x,abs(f(x)-g(x)),'b','LineWidth',4)

% Demonstrate Pythagorean Theorem
l=-1; u=1;
f = @(x) 2*x.^2-1;
g = @(x) 4*x.^3-3*x;
ifsq      = integral(@(x)f(x).^2,l,u);
igsq      = integral(@(x)g(x).^2,l,u);
ifplusgsq = integral(@(x)(f(x)+g(x)).^2,l,u);

% SAVE FIGURES
printfigures(mfilename)