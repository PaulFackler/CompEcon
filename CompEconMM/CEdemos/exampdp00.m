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
h = @(s) s+gamma*(smax-s);


%% SOLUTION

% Code the approximant and the residual
vhat = @(c,s) c(1)+c(2)*s;
vhat1 = @(c,s) price*s-kappa+delta*vhat(c,h(0));
vhat0 = @(c,s) delta*vhat(c,h(s));
resid = @(c,s) vhat(c,s) ...
- max(vhat0(c,s),vhat1(c,s));



% Solve collocation equation
snodes = [0.2;0.4];
c = zeros(2,1);
c = broyden(resid,c,snodes);

% Compute critical biomass
scrit = broyden(@(s) vhat0(c,s)-vhat1(c,s),0)

% delta
% delta*h(snodes(1))
% delta
% delta*h(0)
% price*snodes(1)-kappa

delta
delta*h(snodes(2))
delta
delta*h(0)
price*snodes(2)-kappa

c