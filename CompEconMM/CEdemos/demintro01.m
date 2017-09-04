%% DEMINTRO01 Inverse Demand Problem
%
% Plots demand function q(p)=0.5*p^-0.2+0.5*p^-0.5 and its inverse, and
% computes price that will clear the market of a quantity 2,

% Preliminary tasks
demosetup(mfilename)

% Target quantity and plotting limits
qstar = 2;
pmin  = 0.02;
pmax  = 0.40;

% Compute price that clears market of q=2 using Newton's method
p = 0.25;
for it=1:100
  f = 0.5*p^-0.2 + 0.5*p^-0.5 - qstar;
  d = -0.01*p^-1.2 - 0.25*p^-1.5;
  s = -f/d;
  p = p + s;
  fprintf('iteration %3i  price %8.4f\n',[it p])
  if norm(s)<1.e-8, break, end
end
pstar = p;

% Generate demand function
n = 100;
p = nodeunif(n,pmin,pmax);
q = 0.5*p.^-0.2 + 0.5*p.^-0.5;

% Graph demand and inverse demand functions
figure

subplot(1,2,1)
plot(p,q)
set(gca,'xtick',[])
set(gca,'ytick',[])
set(gca,'box','off')
axis('square')
title('Demand')
xlabel('$p$')
ylabel('$q$')

subplot(1,2,2)
hold on
plot(q,p)
plotvdash(qstar,pstar)
plothdash(qstar,pstar)
plotbullet(qstar,pstar)
plottext(0,pstar,'$p^*$','right','middle');
plottext(qstar,0,'$2$','center','top',14);
set(gca,'xtick',[])
set(gca,'ytick',[])
set(gca,'box','off')
axis('square')
title('Inverse Demand')
xlabel('$q$')
ylabel('$p$')

% Graph market clearing price equilibrium
figure
hold on
f = 0.5*p.^-0.2 + 0.5*p.^-0.5 - 2;
plot(p,f,'LineWidth',5)
plothdash([],0)
plotvdash(pstar,0)
plotbullet(pstar,0)
plottext(0.07,0.85,'$f(p)$')
plottext(pstar+0.005,[],'$p^*$')
set(gca,'XTick',[])
set(gca,'YTick',[0])
set(gca,'YTickLabel',{'0'})
xlabel('$p$')


%% Save Figures
printfigures(mfilename,1)