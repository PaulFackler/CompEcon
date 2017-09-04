%% DEMDDP04 Binomial American Put Option Model

% Preliminary tasks
demosetup(mfilename)


%% FORMULATION

% Model parameters
alpha = 0.2;                            % relative risk aversion
T      = 0.5;                           % years to expiration
sigma  = 0.2;                           % annual volatility
r      = 0.05;                          % annual interest rate
strike = 2.1;                           % option strike price
p0     = 2.0;                           % current asset price

% Discretization parameters
N     = 100;                            % number of time intervals
tau   = T/N;                            % length of time intervals
delta = exp(-r*tau);                    % discount factor
u     = exp(sigma*sqrt(tau));           % up jump factor
q     = 0.5+sqrt(tau)*(r-(sigma^2)/2)/(2*sigma); % up jump probability

% State space
price = p0*(u.^(-N:N))';                % asset prices
n     = length(price);                  % number of states

% Action space (hold=1, exercise=2)
X = (1:2)';                             % vector of actions
m = length(X);                          % number of actions

% Reward function
f = [zeros(n,1) strike-price];

% State transition probability matrix
P = zeros(m,n,n);
for i=1:n
  P(1,i,min(i+1,n)) = q;
  P(1,i,max(i-1,1)) = 1-q;
end

% Model structure
clear model
model.reward    = f;
model.transprob = P;
model.horizon   = N;
model.discount  = delta;


%% SOLUTION

% Solve Bellman Equation
[v,x] = ddpsolve(model);
   

%% ANALYSIS

% Plot optimal exercise boundary
figure
[i,j]=find(diff(x(1:end-1,:)));
plot(flipud((j-1)*tau),price(i))
title('Put Option Optimal Exercise Boundary');
xlabel('Time to Maturity')
ylabel('Asset Price')

% Plot option premium vs. asset price
figure
hold on
plot(price,v(:,1))
plot([0;strike],[strike;0],'k:')
axis([0 strike*2 -inf inf])
title('Put Option Value')
xlabel('Asset Price')
ylabel('Premium')


%% SAVE FIGURES
printfigures(mfilename)