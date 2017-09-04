%% DEMDDP02 Asset Replacement Model

% Preliminary tasks
demosetup(mfilename)


%% FORMULATION

% Model parameters
maxage  = 5;                            % maximum machine age
repcost = 75;                           % replacement cost
delta   = 0.9;                          % discount factor

% State space
S = (1:maxage)';                        % machine age
n = length(S);                          % number of states

% Action space (keep=1, replace=2)
X = (1:2)';                             % vector of actions
m = length(X);                          % number of actions

% Reward function
f = [50-2.5*S-2.5*S.^2 (50-repcost)*ones(n,1)];
f(end,1) = -inf;

% State transition function
g = zeros(n,m);
for i=1:n
  g(i,1) = min(i+1,n);                  % keep
  g(i,2) = 1;                           % replace
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

% Plot value function
figure
plot(S,v)
title('Value Function')
xlabel('Age of Machine')
ylabel('Value of Operation')


%% SIMULATION

% Simulation parameters
nyrs = 12;                              % number of years simulated

% Initialize simulation
sinit = min(S);

% Simulate model
spath = ddpsimul(pstar,sinit,nyrs);

% Plot simulated state path
figure
plot(0:nyrs,S(spath))
title('Simulated Machine Age')
xlabel('Year')
ylabel('Age of Machine')


%% SAVE FIGURES
printfigures(mfilename)