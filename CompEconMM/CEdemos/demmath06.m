%% DEMMATH06 Operations with Markov Chain
%
% Illustrates computation with Markov Chains. 

% Preliminary tasks
demosetup(mfilename)


%% FORMULATION
    
% Model Parameters
gamma  = 0.05;        % aggregate unemployment rate
eta    = 2.0;         % expected duration of unemployment
delta  = 0.9;         % discount factor
y      = [0.5;1.0];   % income per employment state

% Employment Transition Probabilities 
q      = zeros(2,2);
q(1,1) = 1-1/eta;
q(2,1) = gamma*(1-q(1,1))/(1-gamma);
q(1,2) = 1-q(1,1);
q(2,2) = 1-q(2,1);


%% ANALYTIC DERIVATIONS

% Employment State Stationary Distribution and Expected Visit Durations
p   = markov(q);
eta = 1./(1-diag(q));

% Present Value of Expected Income per State
e = (eye(2,2)-delta*q)\y;

% Analytic Ergodic Income Mean, Standard Deviation, and Autocorrelation
Ey = p'*y;
Sy = sqrt(p'*(y.^2)-Ey^2);
Ay = (((p.*y)'*q*y)-Ey^2)/(p'*(y.^2)-Ey^2);


%% SIMULATION

% Number of Replications
n = 100000;

% Perform Simulations
iinit = 1;
isim = [iinit; zeros(n-1,1)];
for j=2:n
  isim(j) = markovsim(isim(j-1),q);
end
ysim = y(isim);

% Visit probabilities
psim = [sum(isim==1)/n sum(isim==2)/n];

% Average length of runs
run = zeros(n,2);
run(1,iinit) = 1;
for j=2:n
  if isim(j)==1
    if isim(j-1)== 1
      run(j,1) = run(j-1,1)+1;
    end
  end
  if isim(j)==2
    if isim(j-1)== 2
      run(j,2) = run(j-1,2)+1;
    end
  end
end
rcount = sum(~run==0);
etasim = sum(run)./rcount;

% Simulated ergodic income mean, standard deviation, and autocorrelation
Eysim = mean(ysim);
Sysim = std(ysim);
Aysim = (sum(ysim(2:n).*ysim(1:n-1))/(n-1)-Eysim^2)/(Sysim^2);


%% COMPARISON OF ANALYTIC AND SIMULATED PROBABILITIES AND MOMENTS

% Print Stationary Distribution and Expected State Durations
fprintf('\n\n')
fprintf('Stationary Distribution and Expected State Durations\n')
fprintf('                   Analytic                 Simulated\n')
fprintf('           Unemployed   Employed    Unemployed   Employed\n')
fprintf('Probability   %7.3f  %9.3f       %7.3f  %9.3f\n',[p(:);psim(:)]')
fprintf('Duration      %7.3f  %9.3f       %7.3f  %9.3f\n',[eta(:);etasim(:)]')

% Print Moments of Ergodic Income Distribution
fprintf('\n\n')
fprintf('Moments of Ergodic Income Distribution\n')
fprintf('                   Analytic Simulated\n')
fprintf('Mean                %7.3f   %7.3f\n',Ey,Eysim)
fprintf('Standard Deviation  %7.3f   %7.3f\n',Sy,Sysim)
fprintf('Autocorrelation     %7.3f   %7.3f\n',Ay,Aysim)