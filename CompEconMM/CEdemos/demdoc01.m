%% DEMDOC01 Deterministic Optimal Consumption-Investment Model
%
% Utility maximizing agent must decide how much to consume and how much to
% hold in a riskless asset.
%
% State
%     w       stock of wealth
% Control
%     q       consumption rate
% Parameters
%     theta   relative risk aversion
%     r       continuous rate of return on asset
%     rho     continuous discount rate

% Preliminary tasks
demosetup(mfilename)


% Initial state and time horizon
winit = 1;        % initial capital stock
T     = 50;       % time horizon


%% SOLUTION & SIMULATION r>rho

% Model parameters
theta = 2.0;                            % relative risk aversion
r     = 0.08;                           % continuous rate of return on asset
rho   = 0.05;                           % continuous discount rate

% V'>0 iff V''<0 iff sign>0 where
sign = rho-r*(1-theta);
if sign<0
  disp('Invalid Parameters')
end

% Solve ODE
g = @(w) ((r-rho)/theta)*w;
[wsim1,~] = oderk4(g,winit,T);


%% SOLUTION & SIMULATION r<rho

% Model Parameters
theta = 2.0;                            % relative risk aversion
r     = 0.05;                           % continuous rate of return on asset
rho   = 0.08;                           % continuous discount rate

% Assume theta>0. Then V'>0 iff V''<0 iff sign>0 where
sign = rho-r*(1-theta);
if sign<0
  disp('Invalid Parameters')
end

% Solve ODE
g = @(w) ((r-rho)/theta)*w;
[wsim2,t] = oderk4(g,winit,T);


%% PLOT SOLUTIONS

% Plot optimal wealth path
figure
plot(t,wsim1,t,wsim2)
legend('$r>\rho$','$\rho>r$')
title('Simulated Wealth')
xlabel('Time')
ylabel('Wealth')


%% SAVE FIGURES
printfigures(mfilename)