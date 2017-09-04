%% DEMDP01A Timber Harvesting Model - Simple Linear Approximation Using Broyden
%
% Profit maximizing a commercial tree stand owner must decide when to
% clear-cut and replant.  This program uses a simple linear approximation
% for the value function and solves the collocation equation direclly using
% Broyden's method.
%
% States
%     s       stand biomass
% Actions
%     j       clear cut/replant (2), not clear cut (1)
% Parameters
%     price   unit price of biomass
%     kappa   clearcut-replant cost
%     smax    stand carrying capacity
%     gamma   biomass growth parameter
%     delta   discount factor

% Preliminary tasks
demosetup(mfilename)


%% FORMULATION
  
% Model parameters
price = 1.0;                           	% unit price of biomass
kappa = 0.2;                            % clearcut-replant cost
smax  = 0.5;                            % stand carrying capacity
gamma = 0.1;                            % biomass growth parameter
delta = 0.9;                            % discount factor

% Growth function
h = @(s) s+gamma*(smax-s);

% Value function approximants and residual
vhat  = @(c,s) c(1)+c(2)*s;
vhat0 = @(c,s)               delta*vhat(c,h(s));
vhat1 = @(c,s) price*s-kappa+delta*vhat(c,h(0));
resid = @(c,s) vhat(c,s) - max(vhat0(c,s),vhat1(c,s));


%% SOLUTION

% Solve collocation equation with two collocation nodes
snodes = [0.2;0.4]; 
c = zeros(2,1);
c = broyden(resid,c,snodes);

% Critical biomass and value
scrit = broyden(@(s) vhat0(c,s)-vhat1(c,s),0);
vcrit = vhat0(c,scrit);
fprintf('Critical Biomass           %9.3f\n',scrit) 


%% ANALYSIS

% Refined state grid
s = nodeunif(1000,0,smax);

% Plot action-contingent value functions
figure
hold on
plot(s,[vhat0(c,s) vhat1(c,s)])
legend('Grow','Clear-Cut','Location','NW')
title('Action-Contingent Value Functions')
xlabel('Biomass')
ylabel('Value of Stand')

% ... plot critical biomass
plotvdash(scrit,vcrit)
plotbullet(scrit,vcrit)
plottext(scrit+0.01,-0.2,'$s^*$')

% Plot residual
figure
hold on
plot(s,100*resid(c,s)./vhat(c,s))
plot(s,0*s,'k-','LineWidth',1)
title('Bellman Equation Residual')
xlabel('Biomass')
ylabel('Percent Residual')

% ... plot collocation nodes
plottext(snodes(1)+0.01,[],'$s_1$')
plottext(snodes(2)+0.01,[],'$s_2$')
plotvdash(snodes(1),0)
plotvdash(snodes(2),0)


%% SIMULATION

% Simulation parameters
nper = 31;                              % number of periods simulated
time = 0:nper-1;                        % periods simulated

% Initialize simulation
s = 0;                                  % initial biomass

% Simulate model
ssim = zeros(nper,1);
for ip=1:nper
  ssim(ip) = s;
  if s<scrit,
    s = h(s);
  else
    s = 0;
  end
end

% Compute mean annual harvest and optimal rotation cycle
s = 0;
for n=1:100
  if s>scrit, break, end
  s = h(s);
end
fprintf('Mean Annual Harvest        %9.3f\n',s/n) 
fprintf('Rotation Cycle in Years    %9i\n\n',n) 

% Plot simulated state path
figure
plot(time,ssim)
title('Simulated Biomass')
xlabel('Period')
ylabel('Biomass')


%% SAVE FIGURES
printfigures(mfilename)