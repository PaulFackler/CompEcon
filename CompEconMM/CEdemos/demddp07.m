%% DEMDDP07 Renewable Resource Model

% Preliminary tasks
demosetup(mfilename)


%% FORMULATION

% Model parameters
delta =  0.9;                           % discount factor
alpha =  4.0;                           % growth function parameter
beta  =  1.0;                           % growth function parameter
gamma =  0.5;                           % demand function parameter
cost  =  0.2;                           % unit cost of harvest

% State space
smin = 0;                               % minimum state
smax = 8;                               % maximum state
n    = 401;                             % number of states
S    = nodeunif(n,smin,smax);           % vector of states

% Action space
xmin = 0;                               % minimum action
xmax = 6;                               % maximum action
m    = 301;                             % number of actions
X    = nodeunif(m,xmin,xmax);           % vector of actions

% Reward function
f = zeros(n,m);
for k=1:m
  f(:,k) = (X(k).^(1-gamma))/(1-gamma)-cost*X(k);
  f(S<X(k),k) = -inf;
end

% State transition function
g = zeros(n,m);
for i=1:n
  for k=1:m
    snext  = alpha*(S(i)-X(k)) - 0.5*beta*(S(i)-X(k)).^2;
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

% Plot optimal policy
figure
plot(S,X(x))
title('Optimal Harvest Policy')
xlabel('Resource Stock')
ylabel('Quantity Harvested')

% Plot value function
figure
plot(S,v)
title('Value Function')
xlabel('Resource Stock')
ylabel('Social Welfare')


%% SIMULATION

% Simulation parameters
nyrs = 20;                              % number of years simulated

% Simulate model
spath = ddpsimul(pstar,n,nyrs);

% Plot simulated state path
figure
plot(0:nyrs,S(spath))
title('Simulated Resource Stock')
xlabel('Year')
ylabel('Resource Stock')

% Plot optimal transition function
figure
[ii,jj]=find(pstar);
plot(S(ii),S(jj),S,S,'--')
title('Optimal State Transition Function')
xlabel('$S_t$')
ylabel('$S_{t+1}$')


%% SAVE FIGURES
printfigures(mfilename)