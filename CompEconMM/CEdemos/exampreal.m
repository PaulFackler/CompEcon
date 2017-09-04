%% Incomplete script file for real analysis examples in lecture notes.

close all


%% Univariate 1st and 2nd order Taylor series approximations

% x = nodeunif(100,-1,1);
% y = (x+1).*exp(2*x);


%% Bivariate 1st and 2nd order Taylor series approximations

% [x1,x2] = ndgrid(0:0.1:2,-1:0.01:1);
% f  = (x1.^2).*exp(-x2);
% f1 = 2*x1 - x2 - 1;
% f2 = x1.^2 - 2*x1.*x2 + 0.5*x2.^2 + x2;

%% Inverse function theorem

% f  = @(p)  0.50*p.^-0.2 + 0.500*p.^-0.5 + 1;
% f1 = @(p) -0.10*p.^-1.2 - 0.250*p.^-1.5;
% f2 = @(p)  0.12*p.^-2.2 + 0.375*p.^-2.5;
% p = nodeunif(100,0.4,3.0);
% q = f(p);


%% Inner product on function space

% l=-1; u=1;
% f = @(x) 2*x.^2-1;
% g = @(x) 4*x.^3-3*x;


%% Norms on Euclidian space

% x = [3;-4];


%% Norm on function space

% l=0; u=2;
% f = @(x) x.^2-1;


%% Metrics on Euclidian space

% x = [7;0];
% y = [4;4];


%% Norm on function space

% l=0; u=2;
% f = @(x) x.^3+x.^2+1;
% g = @(x) x.^3+2;


%% Demonstrate Pythagorean Theorem

% l=-1; u=1;
% f = @(x) 2*x.^2-1;
% g = @(x) 4*x.^3-3*x