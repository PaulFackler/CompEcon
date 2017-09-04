%% DEMDP05 American Put Option Pricing Model
%
% Compute the critical exercise price for an American put option in terms 
% of time to expiration.
%
% States
%     p       underlying asset price
% Actions
%     j       exercize (2) or do not exercize (1)
% Parameters
%     K       option strike price
%     N       number of periods to expiration
%     mu      mean of log price innovation
%     sigma   asset price volatility
%     delta   discount factor

% Preliminary tasks
demosetup(mfilename)


%% FORMULATION

% Continuous time discretization
% sigma = 0.2;                            % annual volatility
% T     = 0.5;                            % years to expiration
% K     = 1.0;                            % option strike price
% r     = 0.1;                            % annual interest rate
% N     = 300;                            % number of time intervals
% dt    = T/N;                            % length of time interval
% delta = exp(-r*dt);                     % per period discount factor
% mu    = dt*(r-sigma^2/2);               % mean of log price innovation
  
% Model parameters
K     = 1.0;                            % option strike price
N     = 300;                            % number of periods to expiration
mu    = 0.0001;                         % mean of log price innovation
sigma = 0.0080;                         % asset price volatility
delta = 0.9998;                         % discount factor
  
% Continuous state shock distribution
m = 15;                                 % number of shocks
[e,w] = qnwnorm(m,mu,sigma^2);          % shocks and probabilities  
  
% Approximation structure
n    = 500;                                     % number of collocation nodes
pmin  = -1;                                     % minimum state
pmax  =  1;                                     % maximum state
[basis,Phi,p] = fundefn('spli',n,pmin,pmax);    % basis functions, interpolation matrix, collocaton nodes


%% SOLUTION
  
% Intialize value function approximant coefficients
c = zeros(n,N+1);

% Solve collocation equation by backward recursion
for t=N:-1:1
  v = zeros(n,1);
  for k=1:m
    pnext = p + e(k);
    v = v + w(k)*max(K-exp(pnext),delta*funeval(c(:,t+1),basis,pnext));
  end
  c(:,t) = Phi\v;
end

% Critical exercise prices
optset('broyden','showiters',0)
pcrit = zeros(N+1,1);
f = @(p,K,delta,c,basis) K-exp(p)-delta*funeval(c,basis,p);
for t=1:N
  pcrit(t) = broyden(f,0,K,delta,c(:,t),basis);
end

% Critical exercise price 300 periods to expiration
fprintf('Critical Exercise Price 300 Periods to Expiration\n')
fprintf('   Critical Price  = %5.2f\n\n' ,exp(pcrit(1)))

% Plot critical exercise prices
figure
time = (N:-1:0);
plot(time,exp(pcrit))
title('American Put Option Optimal Exercise Boundary')
xlabel('Periods Remaining Until Expiration')
ylabel('Exercise Price')


%% SAVE FIGURES
printfigures(mfilename)