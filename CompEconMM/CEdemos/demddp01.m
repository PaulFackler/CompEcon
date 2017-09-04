%% DEMDDP01 Mine Management Model

% Preliminary tasks
demosetup(mfilename)


%% FORMULATION

% Model parameters
price = 1;                              % price of ore
sbar  = 100;                            % initial ore stock
delta = 0.9;                            % discount factor

% State space
S = (0:sbar)';                          % vector of states
n = length(S);                          % number of states

% Action space
X = (0:sbar)';                          % vector of actions
m = length(X);                          % number of actions

% Reward function
f = zeros(n,m);
for i=1:n
  for k=1:m
    if X(k)<=S(i)
      f(i,k) = price*X(k)-(X(k)^2)./(1+S(i));
    else
      f(i,k) = -inf;
    end
  end
end

% State transition function
g = zeros(n,m);
for i=1:n
  for k=1:m
    snext = S(i)-X(k);
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
title('Optimal Extraction Policy')
xlabel('Ore Stock')
ylabel('Quantity Extracted')

% Plot Value Function
figure
plot(S,v)
title('Value Function')
xlabel('Ore Stock')
ylabel('Value of Mine')


%% SIMULATION

% Simulation parameters
nyrs = 15;                              % number of years simulated

% Initialize simulation
sinit = max(S);

% Simulate model
spath = ddpsimul(pstar,sinit,nyrs);

% Plot simulated state path
figure
plot(0:nyrs,S(spath))
title('Simulated Ore Stock Levels')
xlabel('Year')
ylabel('Ore Stock')


%% SAVE FIGURES
printfigures(mfilename)