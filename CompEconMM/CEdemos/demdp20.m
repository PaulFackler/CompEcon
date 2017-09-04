%% DEMDP20 Lifecycle Consumption-Saving Model
%
% In the lifecycle consumption-saving model, an agent works from period 1
% to period T, earning a known income y(t) in each period t, a portion g of
% which he is required to pay into a pension fund. The agent retires at the
% beginning of period T+1, converting the amount accumulated in his pension
% fund plus any other liquid assets (which may be negative) into an annuity
% that provides him with fixed consumption in perpetuity.  The agent begins
% period 1 with no liquid assets, and his liquid assets and pension grow at
% a per-period rate r. The agent may save an unlimited amount, and may
% borrow up to a fixed proportion k of his current income. In this code,
% the Euler conditions are solved as a nonlinear complementarity problem
% using ncpsolve.
%
% Endogenous Variable
%     x       net ending liquid assets (x<0 implies borrowing)
% Parameters
%     g       pension contribution rate (proportion of income)
%     r       interest rate earned by assets and paid on debt
%     alpha   agent's constant relative risk aversion
%     T       number of working periods
%     tmax    period of maximum income
%     ymax    maximum income (income at t=0 is 1)
%     k       borrowing limit as a proportion of income
%     delta   agent's subjective discount factor
% Derived Parameters
%     R       gross interest rate (1+r)
%     y       income per period

function demdp20

% Preliminary tasks
demosetup(mfilename)


%% FORMULATION
  
% Model parameters
g     = 0.00;                           % pension contribution rate (proportion of income)
r     = 0.05;                           % interest rate earned by assets, paid on debt
alpha = 4;                              % agent's relative risk aversion
T     = 40;                             % number of periods until retirement
k     = 1.0;                            % borrowing limit as a proportion of income
delta = 0.90;                           % agent's subjective discount factor
tmax  = 0.75*T;                         % period of maximum income
ymax  = 1.5;                            % maximum income (income at t=0 is 1)

% Derived parameters
t = (1:T)';                             % time periods
R = 1+r;                                % gross interest rate

% Utility function and derivative
u  = @(c) c.^(1-alpha)/(1-alpha);
uder = @(c) c.^(-alpha);

% Income stream
b = 2*(ymax-1)/(tmax-1);
c = b/(tmax-1);
y = 1 + b*(t-1) - 0.5*c*(t-1).^2;

% Plot income stream
figure
plot(t,y)
title('Agent''s Income Stream')
xlabel('Period')
ylabel('Income')

% Value of pension at retirement
G = g*sum(y.*(R.^(T-t+1)));

% Disposable income stream
y = (1-g)*y;


%% SOLUTION - BASE CASE

% Solve Euler conditions as nonlinear complementarity problem
x = 0.2*y;
[x,e] = ncpsolve(@F,-k*y,inf,x,y,T,delta,R,r,G,uder);

% Consumption
c = R*[0;x(1:T-1)] + (1-g)*y - x;
if any(c<0), disp('WARNING: Negative Consumption'), end

% Plot liquid assets & consumption
figure
hold on
plot(t,c,t,x)
plot(t,-k*y,'r:','LineWidth',2)
plothdash([],0)
legend('Consumption','Liquid Assets','Borrowing Limit')
title('Consumption and Liquid Assets')
xlabel('Period')
ylabel('Quantity')

% Compute certainty equivalent consumption and income
U = sum((delta.^(t-1)).*u(c)) + (delta^T/(1-delta))*u(r*(G+R*x(T)));
ccert = ((1-delta)*(1-alpha)*U)^(1/(1-alpha));

% Print output
fprintf('Liquid Assets at Retirement       %5.2f\n',G+R*x(T))
fprintf('Retirement Income and Consumption %5.2f\n',r*(G+R*x(T)))
fprintf('Certainty Equivalent Consumption  %5.2f\n',ccert)
fprintf('Norm of Euler Equation            %9.2e\n',norm(e))

% Plot log marginal utility of consumption
ud = log(uder(R*[0;x(1:T-1)]+y-x));
figure
plot(1:T,ud)
title('Log Marginal Utility of Consumption')
xlabel('Period')


% %% PARAMETRIC ANALYSIS - BORROWING LIMITS
% 
% % Solve for Different Borrowing Limits
% kvec = [0 1 2];
% N = length(kvec);
% 
% % Solve Euler Conditions
% xinit = 0.1*y;
% x = zeros(T,N);
% for i=1:N
%   x(:,i) = ncpsolve(@F,-kvec(i)*y,inf,xinit,y,T,delta,R,r,G,uder);
% end
% s = R*[zeros(1,length(kvec));x(1:T-1,:)];
% 
% % Plot liquid assets
% figure
% hold on
% plot(t,s)
% plothdash([],0)
% legend('$k=0$','$k=1$','$k=2$')
% title('Liquid Assets - Different Borrowing Limits')
% xlabel('Period')
% ylabel('Liquid Assets')


%% PARAMETRIC ANALYSIS - INTEREST RATES

% Solve for Different Interest Rates
rvec = [0.05 0.10 0.15];
Rvec = 1+rvec;
N = length(rvec);

% Solve Euler Conditions
xinit = 0.1*y;
x = zeros(T,N);
for i=1:N
  G = g*sum(y.*(Rvec(i).^(T-t+1)));
  x(:,i) = ncpsolve(@F,-k*y,inf,xinit,y,T,delta,Rvec(i),rvec(i),G,uder);
end
s = R*[zeros(1,length(rvec));x(1:T-1,:)];

% Plot liquid assets
figure
hold on
plot(t,s)
plothdash([],0)
legend('$r=0.05$','$r=0.10$','$r=0.15$')
title('Liquid Assets - Different Interest Rates')
xlabel('Period')
ylabel('Liquid Assets')


% %% PARAMETRIC ANALYSIS - PENSION CONTRIBUTION RATES
% 
% % Solve for Different Pension Contribution Rates
% gvec = [0.0 0.1 0.2];
% N = length(gvec);
% 
% % Solve Euler Conditions
% xinit = 0.1*y;
% x = zeros(T,N);
% for i=1:N
%   G = gvec(i)*sum(y.*(R.^(T-t+1)));
%   x(:,i) = ncpsolve(@F,-k*y,inf,xinit,y,T,delta,R,r,G,uder);
% end
% s = R*[zeros(1,length(gvec));x(1:T-1,:)];
% 
% % Plot liquid assets
% figure
% hold on
% plot(t,s)
% plothdash([],0)
% legend('$g=0.0$','$g=0.1$','$g=0.2$')
% title('Liquid Assets - Different Pension Fund Contribution Rates')
% xlabel('Period')
% ylabel('Liquid Assets')


%% SAVE FIGURES
% printfigures(mfilename)


%% FUNCTION FILE
function [f,J] = F(x,y,T,delta,R,r,G,uder)
ud = [uder(R*[0;x(1:T-1)]+y-x); (r/(1-delta))*uder(r*(G+R*x(T))) ];
t = (1:T)';
f = -ud(t) + delta*R*ud(t+1);
J = [];