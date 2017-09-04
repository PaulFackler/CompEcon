%% DEMSOC02 Stochastic Portfolio Model
%
% Agent must decide the rate at which to consume and how much to invest in
% risky and riskless asset.
%
% State
%     w       wealth
% Control
%     q       consumption rate
%     x       proportion of wealth invested in risky asset
% Parameters
%     r       rate of return on riskless asset
%     R       expected rate of return on risky asset
%     theta   relative risk aversion
%     sigma   volatility or return on risky asset
%     rho     continuous discount rate

% Preliminary tasks
demosetup(mfilename)


%% FORMULATION

% Model parameters
r     = 0.04;                                   % rate of return on riskless asset
R     = 0.08;                                   % expected rate of return on risky asset
theta = 2.0;                                    % relative risk aversion
sigma = 0.05;                                   % volatility or return on risky asset
rho   = 0.05;                                   % continuous discount rate


%% SOLUTION

A = (theta^theta/(1-theta))*(rho-r*(1-theta)-((1-theta)*(r-R)^2)/(2*theta*sigma^2))^(-theta);

% Assume theta>0. Then V'>0 iff V''<0 iff sign>0 where
sign = rho - (1-theta)*(r + ((r-R)^2)/(2*theta*sigma^2));
if sign<0
  disp(['Inviable Parameters'])
end

% Ancillary functions
q = @(w) ((A*(1-theta)).^(-1/theta))*w;
x = (R-r)/(theta*sigma^2);
g = @(w) r*(1-x)*w+R*x*w-q(w);


%% SIMULATION

% Initial state, time horizon, and number of replications
winit = 1;              % initial wealth
T     = 1;              % time horizon
nrep  = 10000;          % number of replications

% Time discretization
N = 1000;               % number of time intervals
h = T/N;                % length of time intervals
t = (0:h:T)';           % time nodes

% Preallocate output array
wsim = zeros(N+1,nrep);

% Simulate model
z = sqrt(h)*randn(N,nrep)*sigma;
w = winit*ones(1,nrep);
for i=0:N
  wsim(i+1,:,:) = w;
  if i<N
    w = w + g(w)*h + x*w.*z(i+1,:);
  end
end

% Plot optimal wealth paths
figure
plot(t,wsim(:,1:3),t,mean(wsim,2),'k')
title('Simulated and Expected Wealth')
xlabel('Time')
ylabel('Wealth')


%% SAVE FIGURES
printfigures(mfilename)