%% ODESPLOT
%
%  Dynamically plots solution of 2-dimensional ODE
%
%  Usage
%    odeplot(x,x1lim,x2lim)
%  Let
%    n = number of observed values of state process
%  Input
%    x        : n.2 observed values of state process
%    x1lim    : 2.1 lower & upper limits of x1
%    x2lim    : 2.1 lower & upper limits of x2

%  Copyright(c) 1997-2015
%   Mario J. Miranda - miranda.4@osu.edu

function odeplot(x,x1lim,x2lim)

n = size(x,1);
m = max(1,floor(n/100));
for j=m+1:n
  if mod(j,m)==0
    plot(x(j-m:j,1),x(j-m:j,2),'r')
    getframe;
    if x(j-m,1)<x1lim(1)||x(j-m,1)>x1lim(2)||x(j-m,2)<x2lim(1)||x(j-m,2)>x2lim(2), break, end
  end
end
plot(x(:,1),x(:,2),'r')