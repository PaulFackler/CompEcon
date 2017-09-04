%% EXAMPDP01  Timber Harvesting Model
%
% Profit maximizing owner of a commercial tree stand must decide when to
% clear-cut the stand.  
%
% States
%     s       stand biomass
% Actions
%     j       clear cut (2) or do not clear cut (1)

% Preliminary tasks
demosetup(mfilename)


%% FORMULATION
  
% Model Parameters
price = 1.0;                           	% price of biomass
kappa = 0.2;                            % clearcut-replant cost
smax  = 0.5;                            % stand carrying capacity
gamma = 0.1;                            % biomass growth parameter
delta = 0.9;                            % discount factor

% Code the growth function


%% SOLUTION

% Code the approximant and the residual

% Solve collocation equation

% Compute critical biomass


%% ANALYSIS

% Compute refined state grid
s = nodeunif(1000,0,smax);

% Plot Action-Contingent Value Functions
figure
hold on
plot(s,vhat0(c,s),s,vhat1(c,s))
legend('Grow','Clear-Cut','Location','N')
title('Action-Contingent Value Functions')
xlabel('Biomass')
ylabel('Value of Stand')

% Compute and Plot Optimal Biomass Harvest Level
plot([scrit scrit],[-0.2 vcrit],'k--')
plotbullet(scrit,-0.2)
plotbullet(scrit,vcrit)
plottext(scrit,-0.2,'$s^*$')
fprintf('Optimal Biomass Harvest Level = %9.4f\n',scrit) 

% Plot Value Function Residual
figure
hold on
plot(s,100*resid(c,s)./vhat(c,s))
plot(s,0*s,'k--')
title('Bellman Equation Residual')
xlabel('Biomass')
ylabel('Percent Residual')
 
% Compute mean annual harvest and optimal rotation cycle
s = 0;
for n=1:100
  if s>scrit, break, end
  s = h(s);
end
fprintf('Mean Annual Harvest        %9.3f\n',s/n) 
fprintf('Rotation Cycle in Years    %9i\n\n',n) 