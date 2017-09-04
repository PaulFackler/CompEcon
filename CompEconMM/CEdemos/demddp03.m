%% DEMDDP03 Asset Replacement Model with Maintenance

% Preliminary tasks
demosetup(mfilename)


%% FORMULATION

% Model parameters
maxage  = 5;                            % maximum asset age
repcost = 75;                           % replacement cost
mancost = 0;                            % maintenance cost
delta   = 0.9;                          % discount factor

% State space
s1 = (1:maxage)';                       % asset age
s2 = (0:maxage-1)';                     % servicings
S  = gridmake(s1,s2);                   % combined state grid
n  = size(S,1);                         % total number of states

% Action space (no action=1, service=2, replace=3)
X = (1:3)';                             % vector of actions
m = length(X);                          % number of actions

% Reward function
f = zeros(n,m);
q = (50-2.5*S(:,1)-2.5*S(:,1).^2);
f(:,1) = q.*min(1,1-(S(:,1)-S(:,2))/maxage);
f(:,2) = q.*min(1,1-(S(:,1)-S(:,2)-1)/maxage) - mancost;
f(:,2) = (50 - repcost)*ones(n,1);

% State transition function
g = zeros(n,m);
for i=1:n
  g(i,1) = getindex([S(i,1)+1 S(i,2)],S);
  g(i,2) = getindex([S(i,1)+1 S(i,2)+1],S);
  g(i,3) = getindex([1 0],S);
end

% Model structure
clear model
model.reward     = f;
model.transfunc  = g;
model.discount   = delta;


%% SOLUTION

% Solve Bellman equation
[v,x,pstar] = ddpsolve(model);
   

%% SIMULATION

% Simulation parameters
nyrs = 12;                              % number of years simulated

% Initialize simulation
sinit = 1;

% Simulate Model
[spath,xpath] = ddpsimul(pstar,sinit,nyrs,x);

% Plot simulated state path (Age)
figure
plot(0:nyrs,S(spath,1))
title('Simulated Asset Age')
xlabel('Year')
ylabel('Age of Asset')
xlim([0 12])

% Plot simulated state path (Servicings)
figure
plot(0:nyrs,S(spath,2))
title('Simulated Number of Servicings')
xlabel('Year')
ylabel('Number of Servicings')
xlim([0 12])
ylim([0 2.25])


%% SAVE FIGURES
printfigures(mfilename)