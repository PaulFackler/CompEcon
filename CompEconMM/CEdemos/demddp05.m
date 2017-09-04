%% DEMDDP05 Water Management Model

% Preliminary tasks
demosetup(mfilename)


%% FORMULATION

% Model parameters
alpha1 = 14;                            % producer benefit function parameter
beta1  = 0.8;                           % producer benefit function parameter
alpha2 = 10;                            % recreational user benefit function parameter
beta2  = 0.4;                           % recreational user benefit function parameter
maxcap = 30;                            % reservoir capacity
r = [0 1 2 3 4];                        % rain levels
p = [0.1 0.2 0.4 0.2 0.1];              % rain probabilities
delta  = 0.9;                           % discount factor

% State space
S = (0:0.1:maxcap)';                   % vector of states
n = length(S);                          % number of states

% Action space
X = (0:0.1:maxcap)';                   % vector of actions
m = length(X);                          % number of actions

% Reward function
f = zeros(n,m);
for i=1:n
  for k=1:m
    if k>i
      f(i,k) = -inf;
    else
      f(i,k) = alpha1*X(k).^beta1+alpha2*(S(i)-X(k)).^beta2;
    end
  end
end

% State transition probability matrix
P = zeros(m,n,n);
for k=1:m
  for i=1:n
    for j=1:length(r)
      snext = min(S(i)-X(k)+r(j),maxcap);
      inext = getindex(snext,S);
      P(k,i,inext) = P(k,i,inext) + p(j);
    end
  end
end

% Model structure
clear model
model.reward     = f;
model.transprob  = P;
model.discount   = delta;


%% SOLUTION

% Solve Bellman equation
[v,x,pstar] = ddpsolve(model);
   

%% ANALYSIS

% Plot optimal policy
figure
plot(S,X(x))
title('Optimal Irrigation Policy')
xlabel('Reservoir Level')
ylabel('Quantity Released for Irrigation')

% Plot value function
figure
plot(S,v)
title('Value Function')
xlabel('Reservoir Level')
ylabel('Social Welfare')


%% SIMULATION

% Simulation parameters
nyrs = 30;                              % number of years simulated
nrep = 10000;                           % number of replications

% Initialize simulation
sinit = ones(nrep,1);

% Simulate model
spath = ddpsimul(pstar,sinit,nyrs);

% Plot simulated state path
figure
plot(0:nyrs,mean(S(spath)))
title('Simulated Reservoir Level')
xlabel('Year')
ylabel('Reservoir Level')

% Compute steady-state distribution of reservoir level
P = markov(pstar);
figure
bar(S,P,10) 
title('Steady State Distribution')
xlabel('Reservoir Level')
ylabel('Probability')

% steady-state reservoir level
avgstock = P'*S;
fprintf('\nSteady-state Reservoir Level    %8.2f\n',avgstock)


%% SAVE FIGURES
printfigures(mfilename)