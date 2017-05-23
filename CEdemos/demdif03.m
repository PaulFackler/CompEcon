% DEMDIF03 Commercial Fishery Model (from V.L. Smith) 
function demdif03
close all

disp(' ')
disp('DEMDIF03 Commercial Fishery Model (from V.L. Smith)')

beta=2.75; f=0.025; del=10;             % Parameter values 

% Compute the trajectories
maxt=30;                                % ending time period   
n=100;                                  % number of time steps 
svals=0:0.1:1;
kvals=0:0.2:2;
x0=phase(svals,kvals);                  % initial conditions 
x0=x0(:,find(x0(1,:)>0));               % eliminate values with S=0
t=[0:0.025:4.95 5:(maxt-5)/n:maxt]';    % time values
[t,x]=rk4('fdif03',t,x0,[],beta,f,del); % call the ODE solver

% Compute the isoclines
s=0.01:0.01:1;
k1=(1-s)./(1-beta*s.*(1-s));            % S zero-isocline 
k2=(1./sqrt(2*f*s)-1./s)./beta;         % K zero-isocline 

% Plot the phase diagram
figure(1)
plot(squeeze(x(:,1,:)),squeeze(x(:,2,:)),'-')
hold on
plot(s,k1,':');
plot(s,k2,':');
hold off
xlabel('S');
ylabel('K');
title('Phase Diagram for Commercial Fishery Example');
axis([0 1 0 2]);
h=text(.005,.3,'K''<0');
set(h,'fontsize',7)
h=text(.06,.1,'K''>0');
set(h,'fontsize',7)
h=text(.825,.15,'S''>0');
set(h,'fontsize',7)
h=text(.825,.35,'S''<0');
set(h,'fontsize',7)


h=text(0.085,1.125,'A');
set(h,'fontsize',8)
h=text(.325,1.8,'B');
set(h,'fontsize',8)
h=text(0.5,1.5,'C');
set(h,'fontsize',8)

prtfigs(mfilename,'Phase Diagram for Commercial Fishery Model',1)



% PHASE   A utility to generate boundary points for 2D phase plot.
% INPUTS: X,Y vectors of values for the x and y axes.
% OUTPUT: P, a 2xK matrix of initial values that can be passed to
%   an ODE solver (e.g., RK4).
function p=phase(x,y);
  miny=min(y);
  maxy=max(y);
  x=x(:)';
  y=y(:)';
  nx=length(x);
  ny=length(y);
  y=y(2:ny-1);
  ny=ny-2;
  p=[x;(miny+zeros(1,nx))];
  if ny>0 p=[p [(max(x)+zeros(1,ny));y]]; end
  p=[p [fliplr(x);(maxy+zeros(1,nx))]];
  if ny>0 p=[p [(min(x)+zeros(1,ny));fliplr(y)]]; end
