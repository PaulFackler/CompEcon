%% DEMDDP11 Stochastic Optimal Growth Model

% Preliminary tasks
demosetup(mfilename)


%% FORMULATION

% Model parameters
delta =  0.9;                           % discount factor
alpha =  0.2;                           % utility parameter
beta  =  0.5;                           % production parameter
gamma =  0.9;                           % capital survival rate

% State space
smin  =  1.0;                           % minimum state
smax  = 10.0;                           % maximum state
n = 200;                                % number of states
S = nodeunif(n,smin,smax);              % vector of states

% Action space
xmin  = 0.65;                           % minimum action
xmax  = 0.99;                           % maximum action
m =  1500;                              % number of actions
X = nodeunif(m,xmin,xmax);              % vector of actions

% Reward function
f = zeros(n,m);
for k=1:m
  f(:,k) = ((S-X(k)*S).^(1-alpha))/(1-alpha);
end

% State transition function
g = zeros(n,m);
for i=1:n
  for k=1:m
    snext = gamma*X(k)*S(i) + (X(k)*S(i)).^beta;
    g(i,k) = getindex(snext,S);
  end
end

% Model structure
clear model
model.reward     = f;
model.transfunc  = g;
model.discount   = delta;


%% SOLUTION

% Solve Bellman equation
[v,x,pstar] = ddpsolve(model);
   

%% ANALYSIS

% Plot optimal investment policy
figure
plot(S,X(x).*S)
title('Optimal Investment Policy')
xlabel('Wealth')
ylabel('Investment')

% Plot optimal consumption policy
figure
plot(S,S-X(x).*S)
title('Optimal Consumption Policy')
xlabel('Wealth')
ylabel('Consumption')

% Plot value function
figure
plot(S,v)
title('Value Function')
xlabel('Wealth')
ylabel('Social Welfare')


%% SIMULATION

% Simulation parameters
nper = 20;                              % number of periods simulated

% Initialize simulation
sinit = smin;                           % initial wealth

% Simulate Model
st = ddpsimul(pstar,sinit,nper);

% Plot simulated state path
figure
plot(0:nper,S(st))
title('Simulated Wealth')
xlabel('Year')
ylabel('Wealth')


%% SAVE FIGURES
printfigures(mfilename)