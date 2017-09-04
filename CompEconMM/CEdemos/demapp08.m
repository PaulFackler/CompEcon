%% DEMAPP08 Linear Spline Approximation

% Preliminary tasks
demosetup(mfilename)

% Function to be approximated
f = @(x) 50-cos(x.^2/8).*(x-pi+.5).^2;

% Basic Figure Setup
xmin = 0;
xmax = 1.5*pi;
off = 0.05;
n = 401;
x = nodeunif(n,xmin,xmax);
y = f(x);
ymin = min(y);
ymax = max(y);
ywid = ymax-ymin;
ymin = ymin-0.5*ywid;
ymax = ymax+0.1*ywid;

% Plot four node approximation
figure
hold on
nnode = 4;
xnode = nodeunif(nnode,xmin,xmax);
ynode = f(xnode);
z = zeros(n,1);
for i=2:nnode
  j = find(x>=xnode(i-1)&x<=xnode(i));
  z(j) = ynode(i-1)+(x(j)-xnode(i-1))*(ynode(i)-ynode(i-1))/(xnode(i)-xnode(i-1));
end
plot(x,y)
plot(x,z,'r','linewidth',2)
xlim([xmin-off xmax+off])
ylim([ymin ymax])
for i=1:nnode
  plotvdash(xnode(i),ynode(i))
  plotbullet(xnode(i),ynode(i))
end
set(gca,'xtick',[],'ytick',[])
plottext(xnode(1),ymin,'$x_0=a$','center','top')
plottext(xnode(2),ymin,'$x_1  $','center','top')
plottext(xnode(3),ymin,'$x_2  $','center','top')
plottext(xnode(4),ymin,'$x_3=b$','center','top')
legend('$f$','$\hat f$','Location','NW')
title('Linear Spline Approximation - Four Nodes')

% Plot seven node approximation
figure
hold on
nnode = 7;
xnode = nodeunif(nnode,xmin,xmax);
ynode = f(xnode);
z = zeros(n,1);
for i=2:nnode
  j = find(x>=xnode(i-1)&x<=xnode(i));
  z(j) = ynode(i-1)+(x(j)-xnode(i-1))*(ynode(i)-ynode(i-1))/(xnode(i)-xnode(i-1));
end
plot(x,y)
plot(x,z,'r','linewidth',2)
xlim([xmin-off xmax+off])
ylim([ymin ymax])
for i=1:nnode
  plotvdash(xnode(i),ynode(i))
  plotbullet(xnode(i),ynode(i))
end
set(gca,'xtick',[],'ytick',[])
plottext(xnode(1),ymin,'$x_0=a$','center','top')
plottext(xnode(2),ymin,'$x_1  $','center','top')
plottext(xnode(3),ymin,'$x_2  $','center','top')
plottext(xnode(4),ymin,'$x_3  $','center','top')
plottext(xnode(5),ymin,'$x_4  $','center','top')
plottext(xnode(6),ymin,'$x_5  $','center','top')
plottext(xnode(7),ymin,'$x_6=b$','center','top')
legend('$f$','$\hat f$','Location','NW')
title('Linear Spline Approximation - Seven Nodes')

%% SAVE FIGURES
printfigures(mfilename)