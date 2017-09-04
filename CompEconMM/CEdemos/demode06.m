%% DEMODE06 Predator-Prey Model
%
%  Solve
%    x1' = -a1*x1         + a2*x1*x2 - h
%    x2' =        + b1*x2 - b2*x1*x2
%  where
%    x1: Predator population
%    x2: Prey population

% Preliminary tasks
demosetup(mfilename)


%% FORMULATION

% Model Parameters
a1 = 0.40;      % natural death rate of predators
a2 = 0.01;      % growth rate of predators per predation
b1 = 1.00;      % natural growth rate of prey
b2 = 0.02;      % death rate of prey per predation
h  = [];        % human hunting of predators

for h=[0 3]     % execute program twice, with hman hunting h = 0 and 3

% Velocity Function
f = @(x) [-a1*x(1,:)+a2*x(1,:).*x(2,:)-h;b1*x(2,:)-b2*x(1,:).*x(2,:)];


%% STEADY-STATE

% Steady State
xst = [b1/b2;(h*b2+a1*b1)/(a2*b1)];
disp('Steady State')
disp(xst)
J = fdjac(f,xst);
disp('Eigenvalues')
disp(eig(J))


%% PHASE DIAGRAM

% Layout
x1lim = [0  80];
x2lim = [0 120];
figlayout(x1lim,x2lim, ...
  'Predator-Prey Model Phase Diagram', ...
  'Predator Population', ...
  'Prey Population')

% Nullclines
x1 = nodeunif(100,x1lim(1),x1lim(2));
x1null = (a1*x1+h)./(a2*x1);
plot(x1,x1null,[b1/b2 b1/b2],x2lim,'LineWidth',2)
legend('Predator Nullcline','Prey Nullcline')

% Velocity Field
odefield(f,x1lim,x2lim)


%% SOLVE ODE USING ODERK4

% Set time discretization
N  = 200;           % number of time nodes

% Solve for given initial value
T  = 30;            % time horizon
x0 = [40;40];       % initial populations
x = oderk4(f,x0,T,N);

% Plot state path
plotbullet(x(1,1),x(1,2),18,'r')
odeplot(x,x1lim,x2lim)
plotbullet(xst(1),xst(2))


%% SOLVE ODE USING ODECOL

% Solve ODE
n = 101;             % number of basis functions
[x,t,res] = odecol(f,x0,T,n);
 
% Plot Solution
figure
plot(t,x)
title('Predator-Prey Model')
xlabel('Time')
ylabel('Population')
legend('Predator','Prey')
set(legend,'Location','NorthWest')

% Plot Residual
figure
hold on
plot(t,res)
plothdash([],0)
title('Predator-Prey Model')
xlabel('Time')
ylabel('Residual')
legend('Predator','Prey')
set(legend,'Location','NorthWest')

end


%% SAVE FIGURES
printfigures(mfilename)