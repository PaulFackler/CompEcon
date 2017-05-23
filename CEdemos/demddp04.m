% DEMDDP04 Binomial American Option Pricing Model
  fprintf('\nDEMDDP04 BINOMIAL AMERICAN OPTION PRICING MODEL\n')
  close all

% Enter model parameters
  T      = 0.5;                 % years to expiration
  sigma  = 0.2;                 % annual volatility
  r      = 0.05;                % annual interest rate
  strike = 2.1;                 % option strike price
  p0     = 2.0;                 % current asset price 

% Enter discretization parameters
  N     = 100;                                     % number of time intervals
  tau   = T/N;                                     % length of time intervals
  delta = exp(-r*tau);                             % discount factor
  u     = exp(sigma*sqrt(tau));                    % up jump factor
  q     = 0.5+sqrt(tau)*(r-(sigma^2)/2)/(2*sigma); % up jump probability
  
  % q=(1/delta-1/u)/(u-1/u); % Alternative value of the up probability

% Construct state space  
  price = p0*(u.^(-N:N))';      % asset prices
  n     = length(price);        % number of states

% Construct reward function (actions hold=1, exercise=2)
  f = [zeros(n,1) strike-price];

% Construct state transition probability matrix
  P = zeros(2,n,n);
  for i=1:n
    P(1,i,min(i+1,n)) = q;
    P(1,i,max(i-1,1)) = 1-q;
  end
  
% Pack model structure
  clear model
  model.reward    = f;
  model.transprob = P;
  model.horizon   = N;
  model.discount  = delta;

% Solve finite-horizon model using backward recursion
  [v,x] = ddpsolve(model);

% Plot option premium vs. asset price
  figure(1); plot(price,v(:,1),'k',[0;strike],[strike;0],'k'); 
  axis([0 strike*2 -inf inf]);
  title('Put Option Value');
  xlabel('Asset Price'); ylabel('Premium');

% Plot optimal exercise boundary
  figure(2); 
  [i,j]=find(diff(x(1:end-1,:)));
  plot(flipud((j-1)*tau),price(i));
  title('Put Option Optimal Exercise Boundary');
  xlabel('Time to Maturity'); ylabel('Asset Price');
    
% Save Plots as EPS Files
%  prtfigs(mfilename,'Solution to the Option Pricing Model',[1 2])