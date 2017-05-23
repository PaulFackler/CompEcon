% DEMQUA02 Compute expectation of function of random normal vector
disp('DEMQUA02 Compute expectation of function of random normal vector')
close all

mu    = [0; 0];                  % enter mean vector
sigma = [1 0.5; 0.5 1];          % enter variance matrix
f = inline('x(:,1).^2 + 2*x(:,1).*x(:,2) - 3*x(:,2).^2');

% Monte Carlo integration
n = 50000;                       % enter number of replicates
x = montnorm(n,mu,sigma);        % generate pseudo-random sequence
yexp = sum(f(x))/n;              % perform monte carlo integration of f
fprintf('Monte Carlo Integration:  %10.3f\n',yexp)

% Gaussian integration
n = [3 3];                       % enter order of approximation
[x,w] = qnwnorm(n,mu,sigma);     % compute Gaussian normal nodes and weights
yexp = w'*f(x);                  % perform Gaussian integration of f
fprintf('Guassian Quadrature:      %10.3f\n',yexp)
