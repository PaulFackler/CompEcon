%% DEMDP16 Linear-Quadratic Model
%
% Solve generic linear-quadratic control model using lqsolve.
%
% States
%     s       generic state of dimension ds=2
% Actions
%     x       generic action of dimension dx=2

% Preliminary tasks
demosetup(mfilename)


%% One-Dimensional Problem

% Model parameters
F0  =  0.0;
Fs  = -1.0;
Fx  = -0.0;
Fss = -1.0;
Fsx =  0.0;
Fxx = -0.1;
G0  =  0.5;
Gs  = -0.2;
Gx  =  0.5;
delta = 0.9; 

% Solve model using lqsolve
[~,~,~,~,~,~,vstar] = lqsolve(F0,Fs,Fx,Fss,Fsx,Fxx,G0,Gs,Gx,delta);

% Alternatively, solve model directly
a = delta*Gx^2;
b = Fxx - delta*(Fxx*Gs^2+Fss*Gx^2-2*Fsx*Gs*Gx);
c = Fsx^2 - Fss*Fxx;
P = (-b-sqrt(b^2-4*a*c))/(2*a);
Hsx = delta*Gs*P*Gx + Fsx;
Hxx = delta*Gx*P*Gx + Fxx;
X = -Hsx/Hxx;
SS = [Fsx Fxx delta*Gx; Fss Fsx delta*Gs-1; Gs-1 Gx  0]\[-Fx;-Fs; -G0];
sstar = SS(1);
xstar = SS(2);
pstar = SS(3);

% Create plot grid
n = 100;                    % number of graphing nodes
smin = -5;                  % minimum state
smax =  5;                  % maximum state
s = nodeunif(n,smin,smax);  % ploting nodes

% Analytic optimal policy and value functions
x = xstar + X*(s-sstar);
v = vstar + pstar*(s-sstar) + 0.5*P*(s-sstar).^2;

% Plot optimal policy
figure
hold on
plot(s,x)
plotbullet(sstar,xstar)
title('Optimal Policy')
xlabel('$s$')
ylabel('$x$')

% Plot value function
figure
hold on
plot(s,v)
plotbullet(sstar,vstar)
title('Value Function')
xlabel('$s$')
ylabel('Value')


%% Higher Dimensional Problem

% Model parameters
ds  = 2;
dx  = 2;
F0  = 3;
Fs  = [1 0];
Fx  = [1 1];
Fss = [-7 -2; -2 -8];
Fsx = [0 0; 0 1];
Fxx = [-2 0; 0 -2];
G0  = [1; 1];
Gs  = [-1 1; 1 0];
Gx  = [-1 -1; 2 3];
delta = 0.95; 

% % Generate Random Model Parameters
% ds  = 2;
% dx  = 2;
% F0 = randn;
% Fs = randn(1,ds);
% Fx = randn(1,dx);
% Fss = randn(ds,ds);
% Fss = -Fss*Fss';
% Fsx = 0.0*randn(ds,dx);
% Fxx = randn(dx,dx);
% Fxx = -Fxx*Fxx';
% G0 = randn(ds,1);
% Gs = randn(ds,ds);
% Gs = -0.1*(Gs*Gs');
% Gx = randn(ds,dx);
% delta = 0.95; 

% Solve model using lqsolve
[X,P,G,sstar,xstar,pstar,vstar] = lqsolve(F0,Fs,Fx,Fss,Fsx,Fxx,G0,Gs,Gx,delta);

% Create plot grid
n = 25 + zeros(1,ds);           % number of graphing nodes per dimension
smin = -1+zeros(1,ds);          % minimum states
smax =  1+zeros(1,ds);          % maximum states
s    = nodeunif(n,smin,smax);   % ploting nodes
s1   = reshape(s(:,1),n);
s2   = reshape(s(:,2),n);

% Derive optimal policy
ns = size(s,1);
x = zeros(ns,dx);
for ix=1:dx
  x(:,ix) = xstar(ix);
  for is=1:ds
    x(:,ix) = x(:,ix) + (s(:,is)-sstar(is))*X(ix,is);
  end
end
x1 = reshape(x(:,1),n);
x2 = reshape(x(:,2),n);

% Derive value function
v = vstar;
for is1=1:ds
  v = v + (s(:,is1)-sstar(is1))*pstar(is1);
  for is2=1:ds
    v = v + 0.5*(s(:,is1)-sstar(is1))*P(is1,is2).*(s(:,is2)-sstar(is2));
  end
end
v  = reshape(v,n);

% Plot Optimal Policy 1
figure
mesh(s1,s2,x1)
title('Optimal Policy')
xlabel('$s_1$')
ylabel('$s_2$')
zlabel('$x_1$')

% Plot Optimal Policy 2
figure
mesh(s1,s2,x2)
title('Optimal Policy')
xlabel('$s_1$')
ylabel('$s_2$')
zlabel('$x_2$')

% Plot Value Function
figure
mesh(s1,s2,v)
title('Value Function')
xlabel('$s_1$')
ylabel('$s_2$')
zlabel('Value')


%% SAVE FIGURES
printfigures(mfilename)