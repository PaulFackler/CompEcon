% exampode02

% Script file for generic ODE example in lecture notes.

% Preliminary tasks
demosetup(mfilename)

% Velocity Function
f = @(x) [x(1,:).^2-2*x(2,:)-5; -1-x(1,:)-x(2,:)];