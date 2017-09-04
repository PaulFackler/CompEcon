%% DEMSOC00 Ito Processes
%
% Simulate geometric Brownian motion

% Preliminary tasks
demosetup(mfilename)

% Model Parameters
T = 1;
n = 365;
h = T/n;
t = (0:h:T)';
mu = 0.1;
sigma = 0.05;

% Simulate
m = 3;
z = randn(n,m);
s = zeros(n+1,m);
s(1,:) = 1;
for i=1:n
  s(i+1,:) = s(i,:) + mu*s(i,:)*h + sigma*s(i,:)*sqrt(h).*z(i,:);
end

% Plot
figure
plot(t,s)
xlabel('$t$')
ylabel('$s(t)$')
title('Simulated Geometric Brownian Motion, $\mu=0.1$, $\sigma=0.05$')

% SAVE FIGURES
printfigures(mfilename)