%% DEMRS03 Dixit Industry Entry-Exit Model
%
% Profit maximizing firm must decide whether to operate or shut down, given
% its current operational state and current profitability, subject to
% transactions costs.

function demrs03

% Preliminary tasks
demosetup(mfilename)


%% FORMULATION

% Model parameters
r     = 0.05;
mu    = 0;
sigma = 0.2;
C     = 1;
E     = 2;
I     = 5;

% Model structure
model.func   = @func;
model.params = {r,mu,sigma,C,E,I};
model.xindex = ...
  [1 0 0 0 0;
  1 2 1 2 0;
  2 1 1 2 0;
  2 0 0 0 1];


%% SOLUTION

% Initial values and approximation size
x = [0; 3; 1; 50];
n = [15 85];

% Solve collocation equation
[cv,basis,x] = rssolve(model,x,n);

x = reshape(x,2,2);
pstar = [x(2,1);x(1,2)];
vstar = [funeval(cv{1},basis{1},pstar(1));
  funeval(cv{2},basis{2},pstar(2))];
dvstar = [funeval(cv{1},basis{1},pstar(1),1);
  funeval(cv{2},basis{2},pstar(2),1)];


%% ANALYSIS

% Value function and derivatives
P = nodeunif(501,0,3);
V = [funeval(cv{1},basis{1},P),...
  funeval(cv{2},basis{2},P)];
V(P>x(2,1),1) = V(P>x(2,1),2)-I;
V(P<x(1,2),2) = V(P<x(1,2),1)-E;
dV = [funeval(cv{1},basis{1},P,1),...
  funeval(cv{2},basis{2},P,1)];
dV(P>x(2,1),1) = dV(P>x(2,1),2);
dV(P<x(1,2),2) = dV(P<x(1,2),1);
d2V = [funeval(cv{1},basis{1},P,2),...
  funeval(cv{2},basis{2},P,2)];
d2V(P>x(2,1),1) = d2V(P>x(2,1),2);
d2V(P<x(1,2),2) = d2V(P<x(1,2),1);
e = r*V-(mu*P(:,[1 1])).*dV-0.5*sigma^2*P(:,[1 1]).^2.*d2V;
e(:,2) = e(:,2)-(P-C);
e(P>x(2,1),1) = nan;
e(P<x(1,2),2) = nan;

% Plot value function
figure
hold on
plot(P,V)
legend('Idle','Active')
title('Value Functions')
xlabel('Profitability')
ylabel('$V$')

%  ... plot switch points
plotvdash(pstar(1),vstar(1))
plotvdash(pstar(2),vstar(2))
plotbullet(pstar,vstar)
plottext(pstar(1)+0.02,[],'$p_h^*$')
plottext(pstar(2)+0.02,[],'$p_l^*$')

% Plot first derivative of value function
figure
hold on
dV(P>pstar(1),1) = NaN;
dV(P<pstar(2),2) = NaN;
plot(P,dV)
legend('Idle','Active')
title('First Derivative of Value Function')
xlabel('Profitability')
ylabel('$V''$')

%  ... plot switch points
plotvdash(pstar(1),dvstar(1))
plotvdash(pstar(2),dvstar(2))
plotbullet(pstar(1),dvstar(1))
plotbullet(pstar(2),dvstar(2))
plottext(pstar(1)+0.02,[],'$p_h^*$')
plottext(pstar(2)+0.02,[],'$p_l^*$')

% Plot residual
figure
hold on
plot(P,e)
plothdash([],0)
legend('Idle','Active')
title('Approximation Residual')
xlabel('Profitability')
ylabel('Residual')

% "Closed-form" solutions
optset('broyden','showiters',0)
optset('broyden','maxit',100)
[~,~,~,~,~,Pl,Ph] = entexgbm([],r,mu,sigma,I,E,C);

fprintf('Switching Points\n') 
fprintf('        "Exact"      Approximate        Error\n') 
fprintf('Ph       %5.3f         %5.3f         %8.1e\n',[Ph pstar(1) Ph-pstar(1)])
fprintf('Pl       %5.3f         %5.3f         %8.1e\n',[Pl pstar(2) Pl-pstar(2)])


%% SAVE FIGURES
printfigures(mfilename)


%% RSSOLVE FUNCTION FILE
%
%   User-supplied function that returns the optimal control, reward, state
%   drift, and state diffusion at an arbitrary number of ns states:
%      out = func(flag,s,x,<params>)

function out = func(flag,s,x,r,mu,sigma,C,E,I)
switch flag
  case 'f'        % reward
    out = (s-C).*(x==2);
  case 'g'        % state drift
    out = mu*s;
  case 'sigma'    % state diffusion
    out = sigma*s;
  case 'rho'      % state contingent discount rate
    out = r+zeros(size(s,1),1);
  case 'reward'   % jump rewards and derivatives
    out = [0 0 0;-I 0 0;-E 0 0;0 0 0];
end


% ENTEXGBM Finds entry and exit triggers for Geometric Brownian Motion
% ref.: Dixit & Pindyck, p. 218
function [residual,beta1,beta2,A1,A2,pl,ph] = ...
  entexgbm(p,rho,mu,sigma,I,E,c,pfactor,zl,zh,beta1,beta2)

if nargin==12   % ******** Computes residuals for ROOT ***********
  pl = exp(p(1,:));  % converted to logs to avoid negative values
  ph = exp(p(2,:));
  A1 = p(3,:);
  A2 = p(4,:);
  
  % Compute residuals from value matching & smooth pasting conditions
  residual = [
    (A1.*ph.^beta1 - A2.*ph.^beta2 - ph*pfactor + zh);
    (A1.*pl.^beta1 - A2.*pl.^beta2 - pl*pfactor + zl);
    (beta1*A1.*ph.^(beta1-1) - beta2*A2.*ph.^(beta2-1) - pfactor);
    (beta1*A1.*pl.^(beta1-1) - beta2*A2.*pl.^(beta2-1) - pfactor)];
  
elseif nargin==7 % ******** Set up problem & call BROYDEN ********
  
  if mu>=rho || sigma<=0 || c<0 || rho<=0 || I<0
    error('Improper parameter values');
  end
  
  pfactor = 1/(rho-mu);
  zh = c/rho+I;
  zl = c/rho-E;
  
  s2 = sigma.^2;
  beta1 = 0.5-mu/s2+sqrt((0.5*s2-mu).^2+2*rho*s2)/s2;
  beta2 = 0.5-mu/s2-sqrt((0.5*s2-mu).^2+2*rho*s2)/s2;
  
  % beta = roots([s2/2 mu-s2/2 -rho]);
  
  if isempty(p), p = [0.1;1;1;1]; end   % initialize P if starting value empty
  
  mfilename
  p = broyden(@entexgbm,log(p),rho,mu,sigma,I,E,c,pfactor,zl,zh,beta1,beta2);
  % compute residuals - should be close to zero
  residual = entexgbm(p,rho,mu,sigma,I,E,c,pfactor,zl,zh,beta1,beta2);
  
  pl = exp(p(1,:));  % converted to logs to avoid negative values
  ph = exp(p(2,:));
  A1 = p(3,:);
  A2 = p(4,:);
  
else
  error('Wrong number of input arguments')
end