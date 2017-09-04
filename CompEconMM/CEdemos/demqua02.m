%% DEMQUA02 Expectation of Function of Random Normal Vector
%
%  Expectation of function of bivariate normal vector tsing Monte Carlo and
%  Gaussian quadrature

% Preliminary tasks
demosetup(mfilename)

% Distribution parameters
mu    = [0;0];                   % mean vector
sigma = [1 0.5; 0.5 1];          % variance matrix

% Function to be integrated
f = @(x) exp(-x(:,1)).*cos(x(:,2).^2); 

% Monte Carlo integration
n = 50000;                       % number of replicates
x = mvnrnd(mu,sigma,n);          % pseudo-random sequence
yexp = sum(f(x))/n;              % monte carlo integration of f
fprintf('Monte Carlo Integration:  %10.3f\n',yexp)

% Gaussian integration
n = [21 21];                     % order of approximation
[x,w] = qnwnorm(n,mu,sigma);     % Gaussian normal nodes and weights
yexp = w'*f(x);                  % Gaussian integration of f
fprintf('Guassian Quadrature:      %10.3f\n',yexp)