%% DEMINTRO02 Rational Expectations Agricultural Market Model
%
% Compute rational expectations equilibrium of agricultural market model
% using function iteration.

% Preliminary tasks
demosetup(mfilename)


%% Generate Yield Distribution

% Yield distribution
sigma = 0.2;                                % Log standard deviation of yield
[y,w] = qnwlogn(25,-0.5*sigma^2,sigma^2);   % Discretizied log-normal distribution, mean 1


%% Compute equilibrium, target price = 1

% Set target price
ptarg = 1;

% Compute equilibrium, iterating on acreage planted
fprintf('\n\n')
a = 1;
for it=1:50
  aold = a;
  Ef = w'*max(1.5-0.5*a*y,ptarg);
  a = 0.5 + 0.5*Ef;
  fprintf('iteration %3i  acreage %8.4f   expected price %8.4f\n',[it a Ef])
  if norm(a-aold)<1.e-8, break, end
end

% % Compute equilibrium, iterating on expected farm price
% fprintf('\n\n')
% Ef = 1;
% for it=1:50
%   Efold = Ef;
%   a = 0.5 + 0.5*Ef;
%   Ef = w'*(max(1.5-0.5*a*y,ptarg));
%   fprintf('iteration %3i  acreage %8.4f   expected price %8.4f\n',[it a Ef])
%   if norm(Ef-Efold)<1.e-8, break, end
% end

% Intermediate output
q = a*y;            % quantity produced in each yield state
p = 1.5-0.5*a*y;    % market price in each yield state
f = max(p,ptarg);   % farm price in each yield state
r = f.*q;           % farm revenue in each yield state
g = (f-p).*q;       % government expenditures in each yield state

% Compute means and standard deviation
x = [p f r g];
xavg = w'*x;
xstd = sqrt(w'*x.^2-xavg.^2);

% Print results 
fprintf('\n\n')
fprintf('Variable                   Expect     Std Dev\n')
fprintf('Market Price             %8.4f   %8.4f\n',xavg(1),xstd(1))
fprintf('Farm Price               %8.4f   %8.4f\n',xavg(2),xstd(2))
fprintf('Farm Revenue             %8.4f   %8.4f\n',xavg(3),xstd(3))
fprintf('Government Expenditures  %8.4f   %8.4f\n',xavg(4),xstd(4))


%% Plot fixed-point

% Generate fixed-point mapping
aeq = a;
a = nodeunif(100,0,2);
g = zeros(100,1);
for i=1:100
  g(i) = 0.5 + 0.5*w'*max(1.5-0.5*a(i)*y,1);
end

% Graph rational expectations equilibrium
figure
hold on
plot(a,g)
plot(a,a,'k:','LineWidth',2)
plotvdash(aeq,aeq)
plothdash(aeq,aeq)
plotbullet(aeq,aeq)
plottext(0.05,1.25,'$g(a)$')
plottext(0.15,0.00,'$45^\circ$',[],[],12)
plottext(-0.15,aeq-0.10,'$a^*$')
plottext(aeq-0.05,-0.20,'$a^*$')
axis('square')
set(gca,'XTick',[0 2],'XTickLabel',{'0' '2'})
set(gca,'YTick',[0 2],'YTickLabel',{'0' '2'})
title('Equilibrium Acreage Planted')
xlabel('Acreage Planted')


%% Sensitivity analysis with respect to target price

% Compute equilibrium as function of target price
nplot = 50;
ptarg = nodeunif(nplot,0,2);
a = 1;
for ip=1:nplot
  for it=1:50
    aold = a;
    a = 0.5 + 0.5*w'*max(1.5-0.5*a*y,ptarg(ip));
    if norm(a-aold)<1.e-10, break, end
  end
  q = a*y;              % quantity produced
  p = 1.5-0.5*a*y;      % market price
  f = max(p,ptarg(ip)); % farm price
  r = f.*q;             % farm revenue
  g = (f-p).*q;         % government expenditures
  x = [p f r g];
  xavg = w'*x;
  xstd = sqrt(w'*x.^2-xavg.^2);
  Ep(ip) = xavg(1);
  Ef(ip) = xavg(2);  
  Er(ip) = xavg(3);
  Eg(ip) = xavg(4);
  Sp(ip) = xstd(1);
  Sf(ip) = xstd(2);  
  Sr(ip) = xstd(3);
  Sg(ip) = xstd(4);
end

% Graph expected prices vs target price
figure
hold on
plot(ptarg,Ep,ptarg,Ef,'LineWidth',5)
plot(ptarg([1 nplot]),Ep([1 1]),'k:','LineWidth',2)
legend('Market Price','Farm Price','Location','NW')
set(gca,'XTick',[0 1 2])
set(gca,'YTick',[0.5 1.0 1.5 2.0])
set(gca,'YTickLabel',{'0.5' '1.0' '1.5' '2.0'})
ylim([0.5 2])
title('Expected Prices vs Target Price')
xlabel('Target Price')
ylabel('Expected Price')

% Graph price standard deviations vs target price
figure
hold on
plot(ptarg,Sp,ptarg,Sf,'LineWidth',5)
plot(ptarg([1 nplot]),Sf([1 1]),'k:','LineWidth',2)
legend('Market Price','Farm Price','Location','NW')
set(gca,'XTick',[0 1 2])
set(gca,'YTick',[0 0.1 0.2])
title('Price Variabilities vs Target Price')
xlabel('Target Price')
ylabel('Standard Deviation of Price')

% Graph expected farm revenue vs target price
figure
hold on
plot(ptarg,Er,'LineWidth',5)
plot(ptarg([1 nplot]),Er([1 1]),'k:','LineWidth',2)
set(gca,'XTick',[0 1 2])
set(gca,'YTick',[0 1 2 3])
set(gca,'YTickLabel',{'0' '1' '2' '3'})
ylim([0.8 3])
title('Expected Farm Revenue vs Target Price')
xlabel('Target Price')
ylabel('Expected Farm Revenue')

% Graph standard deviation of farm revenue vs target price
figure
hold on
plot(ptarg,Sr,'LineWidth',5)
plot(ptarg([1 nplot]),Sr([1 1]),'k:','LineWidth',2)
set(gca,'XTick',[0 1 2])
set(gca,'YTick',[0 0.3 0.6])
set(gca,'YTickLabel',{'0' '0.3' '0.6'})
ylim([0. 0.6])
title('Farm Revenue Variability vs Target Price')
xlabel('Target Price')
ylabel('Standard Deviation of Farm Revenue')

% Graph expected government expenditures vs target price
figure
hold on
plot(ptarg,Eg,'LineWidth',5)
plothdash([],0)
set(gca,'XTick',[0 1 2])
set(gca,'YTick',[0 1 2])
ylim([-0.05 2])
title('Expected Government Expenditures vs Target Price')
xlabel('Target Price')
ylabel('Expected Government Expenditures')


%% Save Figures
printfigures(mfilename)