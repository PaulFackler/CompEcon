%% DEMDP18 Credit With Technology Adoption Model
%
% An infinitely-lived agent subject to production and income uncertatintly
% must decide how much to consume, how much to borrow, and whether to
% invest in a hihgher yielding technolocy. At the beginning of any period,
% the agent possses a disposable wealth w, which may be negative if his
% income is less than his debt obligation. The agent may borrow up to a
% specified limit.  He may also invest in a higher yielding tehcnology.
%
% States
%      w        agent net worth (w<0 indicates debt)
% Actions
%      x        amount borrowed
%      j        technology adoption decision (1=lo tech, 2=hi tech)

function demdp18

% Preliminary tasks
demosetup(mfilename)


%% FORMULATION

% Numerical control parameters
mz = 3;                                 % number of off-farm income shocks
my = 4;                                 % number of on-farm income shocks

% Model parameters
rho    = 0.05;                          % subjective discount rate
r      = 0.05;                          % interest rate on borrowing
R      = 0.08;                          % rate of return on high-yield technology
alpha  = 2.0;                           % relative risk aversion
barz   = 0.50;                          % expected off-farm income
bary   = [1.0;1.1];                     % expected on-farm income, lo tech, hi tech
sigmaz = 0.2;                           % off-farm income volatility
sigmay = 0.2;                           % on-farm income volatility
xlim   = 0.5;                           % borrowing limit

% Derived parameters
delta = 1/(1+rho);                      % discount factor
K = (bary(2)-bary(1))/(1+R);            % cost of high-yield technology

% Income distribution
[e1,p1] = qnwlogn(mz,-0.5*sigmaz^2,sigmaz^2);       % off-farm income shocks
[e2,p2] = qnwlogn(my,-0.5*sigmay^2,sigmay^2);       % on-farm income shocks
e = gridmake(e1,e2);                                % joint income shock distribution shocks
p = kron(p2,p1);                                    % joint income shock distribution probabilities

% Wealth bounds
wmin = min(barz*e(:,1)+bary(1)*e(:,2))-(1+r)*xlim;  % minimum wealth
wmax = max(barz*e(:,1)+bary(2)*e(:,2));             % maximum wealth
if wmin+xlim<0, error('Inadequate Wealth Possible'), end

% Model structure
model.func = @func;                     % model functions
model.params = {alpha barz bary K r xlim};   % function parameters
model.discount = delta;                 % discount factor
model.ds = 1;                           % dimension of continuous state
model.dx = 1;                           % dimension of continuous action
model.ni = 0;                           % number of discrete states
model.nj = 2;                           % number of discrete actions
model.e  = e;                           % continuous state shocks
model.w  = p;                           % continuous state shock probabilities

% Approximation structure
n = 100;                                % number of collocation nodes
basis = fundefn('spli',n,wmin,wmax);  	% basis functions

% Check function files
dpcheck(model,(wmin+wmax)/2,xlim/2,1,1)


%% SOLUTION

% Solve collocation equation
[~,w,v,x,resid] = dpsolve(model,basis);


%% ANALYSIS

% Plot optimal policy
figure
plot(w,x)
legend('Lo Tech','Hi Tech')
title('Action Contingent Optimal Borrowing Policy')
xlabel('Wealth')
ylabel('Amount Borrowed')

% Plot value function
figure
plot(w,v)
legend('Lo Tech','Hi Tech')
title('Action-Contingent Value Functions')
xlabel('Wealth')
ylabel('Lifetime Utility')

% Plot residual
figure
hold on
plot(w,resid)
plothdash([],0)
title('Bellman Equation Residual')
xlabel('Wealth')
ylabel('Residual')


%% SIMULATION

% Simulation parameters
nper = 20;                              % number of periods simulated
nrep = 2000;                            % number of agents replicated

% Initialize simulation
wtmp = wmin*ones(nrep,1);               % all agents have minimal wealth
rng('default')        	

% Pre-allocate history arrays
xtmp = zeros(nrep,1);
wsim = zeros(nrep,nper+1);
xsim = zeros(nrep,nper+1);
jsim = zeros(nrep,nper+1);

% Simulate model
for k=1:nper+1
  wsim(:,k) = wtmp;
  j = discrand(nrep,p);
  v1 = interp1(w,v(:,1),wtmp);
  v2 = interp1(w,v(:,2),wtmp);
  i = find(v2>v1);
  if ~isempty(i)
    xtmp(i) = interp1(w,x(:,2),wtmp(i));
    wtmp(i) = barz*e(j(i),1) + bary(2)*e(j(i),2) - (1+r)*xtmp(i);
    jsim(i,k) = 1;
  end
  i = find(v1>=v2);
  if ~isempty(i)
    xtmp(i) = interp1(w,x(:,1),wtmp(i));
    wtmp(i) = barz*e(j(i),1) + bary(1)*e(j(i),2) - (1+r)*xtmp(i);
  end
  xsim(:,k) = xtmp;
end

% Plot simulated agent wealth
figure
hold on
plot(0:nper,wsim(1:3,:)','LineWidth',2)
plot(0:nper,mean(wsim),'k')
title('Agent''s Simulated and Expected Wealth')
xlabel('Period')
ylabel('Wealth')

% Plot simulated agent borrowing
figure
hold on
plot(0:nper,xsim(1:3,:)','LineWidth',2)
plot(0:nper,mean(xsim),'k')
title('Agent''s Simulated and Expected Borrowing')
xlabel('Period')
ylabel('Amount Borrowed')

% Plot simulated technology adoption
figure
hold on
plot(0:nper,mean(jsim),'k')
title('Technology Adoption')
xlabel('Period')
ylabel('Proportion of Agents Using High Technology')


%% SAVE FIGURES
printfigures(mfilename)


%% DPSOLVE FUNCTION FILE
%
%    User-supplied function called by dpsolve that returns the bound,
%    reward, and continuous state transition function values and
%    derivatives with respect to the continuous action x at an arbitrary
%    number ns of states and actions according to the format
%      [out1,out2,out3] = func(flag,s,x,i,j,e,<params>)
%    where s is ns.ds continuous states, x is ns.dx continuous actions, i
%    is ns.1 or scalar discrete states, j is ns.1 or scalar discrete
%    actions, and e is ns.de continuous state transition shocks.

function [out1,out2,out3] = func(flag,w,x,~,j,e,alpha,barz,bary,K,r,xlim)

cmin = 0.05;
n = length(w);
switch flag
  case 'b'      % bounds
    out1 = max(0,cmin - w + (j-1)*K);
    out2 = xlim+zeros(n,1);
  case 'f'      % reward
    c =  w + x - (j-1)*K;
    out1 = (c.^(1-alpha))/(1-alpha);
    out2 = (c.^(-alpha));
    out3 = -alpha*(c.^(-alpha-1));
  case 'g'      % transition
    out1 = barz*e(:,1) + bary(j)*e(:,2) - (1+r)*x;
    out2 = -(1+r)*ones(n,1);
    out3 = zeros(n,1);
end