%% DEMODE08 Lorentz Strange Attractor
%
%  Solve
%   x1' = -a*x1 +a*x2
%   x2' = -x2 - x1*x3
%   x3' = -b*x3 + x1*x2 - b*c;
%  s.t.
% 	x1(0) = -8
% 	x2(0) =  8
% 	x3(0) = 27

% Preliminary tasks
demosetup(mfilename)



%% FORMULATION

% Model Parameters
a = 10;
b = 8/3;
c = 28;

% Velocity Function
f = @(x) [-a*x(1,:)+a*x(2,:);-x(2,:)-x(1,:).*x(3,:);-b*x(3,:) + x(1,:).*x(2,:)-b*c];


%% SOLVE ODE USING ODERK4

% Solve ODE
N = 750;                    % number of time nodes
T = 15;                     % time horizon
x0 = [-8;8;27];             % initial state
[x,t] = oderk4(f,x0,T,N);

% Plot Solution
figure
plot(t,x)
title('Lorentz Strange Attractor')
xlabel('Time')
ylabel('State')
legend('$x_1$','$x_2$','$x_3$')


%% PHASE DIAGRAM

% Layout
figlayout([-20 20],[-30 30], ...
  'Lorentz Strange Attractor Phase Diagram','$x_1$','$x_3$')

% State Path
for j=2:N
   plot(x(j-1:j,1),x(j-1:j,3),'r','LineWidth',2)
   getframe; 
end
plot(x(:,1),x(:,3),'r','LineWidth',2)

% plotbullet(x1(1,1),x1(1,2),18,'r')
% plotbullet(x2(1,1),x2(1,2),18,'r')
% plotbullet(kst,cst,20)


%% SAVE FIGURES
printfigures(mfilename)
