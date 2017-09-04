%% DEMDP17 Miscellaneous Lecture Note Figures

% Preliminary tasks
demosetup(mfilename)

optset('broyden','showiters',0)


%% FIGURE 1 Timber growth curve

a = 0.1;
sstar = 0.9;
c = 2*(sstar-a)/sstar^2;
b = c*sstar;

% Define timber growth curve
F = @(s) a+b*s-0.5*c*s.^2;

% Figure setup
n = 401;
smin = 0;
smax = 1;
s  = nodeunif(n,smin,smax);

% Plot timber growth curve
figure
hold on
% Figure setup
axis square
xlim ([0 1])
ylim ([0 1])
set(gca,'xtick',[])
set(gca,'ytick',[])
% Plot growth function
plot(s,F(s))
% Plot 45-degree line
plot(s,s,'k:','Linewidth',2)
plottext(0.12,0.01,'$45^\circ$','right','bottom',12)
% Plot fixed points
plotbullet(sstar,sstar)
plotvdash(sstar,sstar)
plottext(0,0.95,'$h(s)$','right','middle')
plottext(1,0-.015,'$s$','center','top')
plottext(0.91,0-.015,'$\bar s$')
plottext(0-0.01,a,'$h(0)$','right','middle',12)


%% FIGURE 2 Asset productivity curves

% Figure setup
n = 401;
amin = 0;
amax = 1;
a  = nodeunif(n,amin,amax);

% Define timber growth curve
f1 = @(a) 0.9-0.4*(exp(a)-1);
f2 = @(a) 0.9-2.7*(a-0.5).^2;

figure
subplot(1,2,1)
hold on
plot(a,f1(a))
xlabel('Age')
ylabel('Output')
title('Declining')
set(gca,'xtick',[])
set(gca,'ytick',[])
xlim ([0 1])
ylim ([0 1])
axis square
subplot(1,2,2)
hold on
plot(a,f2(a))
xlabel('Age')
ylabel('Output')
title('Concave')
set(gca,'xtick',[])
set(gca,'ytick',[])
xlim ([0 1])
ylim ([0 1])
axis square


%% FIGURE 3 consumer surplus

% Initiate figure
figure
hold on

% Define marginal cost curve
f = @(q) 0.2+0.7*exp(-5*q);

% Set plotting parameters
n = 1001;
qmin = 0;
qmax = 1;
pmin = 0;
pmax = 1;
q = nodeunif(n,0,qmax);
p = f(q);

% Plot area under inverse demand curve
qstar = 0.4;
pstar = f(qstar);
kstar = 0.5*pstar;
for i=1:n
  if 0.005<q(i) && q(i)<qstar
    plot([q(i) q(i)],[kstar p(i)-0.02],'y-','linewidth',2)
  end
end

% Plot inverse demand curve
plot(q,p)
xlim([qmin qmax])
ylim([pmin pmax])
set(gca,'xtick',[])
set(gca,'ytick',[])
axis square
% xlabel('Quantity')
% ylabel('Price')
plottext(qmin,pmin,'$0$','right','top')
plottext(qmax,pmin,'$q$','center','top')
plottext(qmin-0.02,pmax,'$p$','right','top')

% Annotate figure
plottext(qmin-0.02,pstar,'$p_t$','right','middle')
plottext(qmin-0.02,kstar,'$k$','right','middle')
plottext(qstar,pmin,'$q_t$','center','top');
plothdash(qstar,pstar)
plothdash(qstar,kstar)
plotvdash(qstar,pstar)
qtmp = 0.4*qstar; ptmp = 0.5*f(qtmp);
plottext(qtmp,ptmp-0.04,'PS','center','middle')
plottext(qtmp,ptmp+0.11,'CS','center','middle')


%% FIGURE 4 Renewable resource model biological growth function

% Gross return
R  = 1.2;

% Define biological growth function
m1 = 0.4;
m2 = 0.6;
s1 = 1.0;
s2 = 1.0;
c  = 0.45;
F1 = @(q) cdf('norm',q,m1,(s1*m1)^2);
F2 = @(q) cdf('norm',q,m2,(s2*m2)^2);
f1 = @(q) pdf('norm',q,m1,(s1*m1)^2);
f2 = @(q) pdf('norm',q,m2,(s2*m2)^2);
c2 = c; c1 = c*(f2(0))/(f1(0));
F  = @(q) c1*(F1(q)-F1(0))-c2*(F2(q)-F2(0));

% Compute fixed points
G  = @(q) F(q)-q;
q1 = broyden(G,m1);
q2 = broyden(G,m2);

% Compute optimal retention, stock, harvest
g  = @(q) c1*f1(q)-c2*f2(q)-R;
optset('broyden','defaults');
rstar = broyden(g,m2);
sstar = F(rstar);

% Figure setup
n = 401;
qmin = 0;
qmax = 1;
q = nodeunif(n,qmin,qmax);

% Plot growth function and steady-states
figure
hold on
% Figure setup
axis square
ylim ([0 1])
set(gca,'xtick',[])
set(gca,'ytick',[])
% Plot growth function
plot(q,F(q))
% Plot 45-degree line
plot(q,q,'k:','Linewidth',2)
plottext(1,0.73,'$g$','right','middle')
plottext(0.08,0.08,'$45^\circ$','right','bottom',12)
% Plot fixed points
plotbullet( 0 ,0)
plotbullet(q1,q1)
plotbullet(q2,q2)
plotvdash(q1,q1)
plotvdash(q2,q2)
yl = ylim;
plottext(q1+0.01,yl(1),'$s^*_1$')
plottext(q2+0.01,yl(1),'$s^*_2$')
xlabel('Stock This Period')
ylabel('Stock Next Period')


%% FIGURE 5 Renewable resource model optimal solution

% Plot growth function and optimal solution
figure
hold on
% Figure setup
axis square
ylim ([0 1])
set(gca,'xtick',[])
set(gca,'ytick',[])
% Plot growth function
plot(q,F(q))
% Plot 45-degree line
plot(q,q,'k:','Linewidth',2)
plottext(0.08,0.08,'$45^\circ$','right','bottom',12)
% Plot slope line
qq = nodeunif(n,rstar-0.2,rstar+0.2);
yy = F(rstar)+R*(qq-rstar);
plot(qq,yy,'k:','linewidth',2)
plot(rstar,sstar,'k*','linewidth',8)
% Plot optimal retention, stock, harvest
fudge = 0.02;
plothdash(rstar,sstar)
plothdash(rstar,rstar)
plotvdash(rstar,rstar)
plotvdash(rstar,sstar)
plottext(1,0.73,'$g$','right','middle')
plottext(0-fudge,rstar,'$s^*-q^*$','right','middle')
plottext(0-fudge,sstar,'$s^*$','right','middle')
plotbullet(rstar,rstar)
plottext(0.65,0.9,'slope=$1+\rho$','right','bottom',10)
plottext(rstar,0,'$s^*-q^*$','center','top')
mid = 0.5*(rstar+sstar);
plottext(0+fudge,mid,'$q^*$','left','middle')
annotation('arrow',[0.23 0.23],[mid-0.04 rstar+0.01],'LineWidth',0.1)
annotation('arrow',[0.23 0.23],[mid+0.01 sstar-0.04],'LineWidth',0.1)


%% FIGURE 6 Cost of Extraction 

% Initiate figure
figure
hold on

% Define marginal cost curve
f = @(s) 0.2+0.7*exp(-5*s);

% Set plotting parameters
n = 401;
smin = 0;
smax = 1;
kmin = 0;
kmax = 1;
s = nodeunif(n,smin,smax);

% Plot marginal cost curve
k = f(s);
plot(s,k)
xlim([smin smax])
ylim([kmin kmax])
set(gca,'xtick',[])
set(gca,'ytick',[])
% xlabel('Ore Stock')
% ylabel('Marginal Cost of Extraction')
plottext(smin,kmin,'$0$','right','top')
plottext(smin,kmax,'$k(s)$','right','top')
plottext(smax,f(smax)-0.2,'$s$','center','top')

% Plot area under marginal cost curve
a = 0.50;
b = 0.75;
for i=1:n
  if a<=s(i) && s(i)<=b
    plot([s(i) s(i)],[kmin+0.005 k(i)-0.01],'y-','linewidth',2)
  end
end
% plot([a a],[kmin f(a)],'k--','linewidth',2)
% plot([b b],[kmin f(b)],'k--','linewidth',2)
plotvdash(a,f(a))
plotvdash(b,f(b))
plottext(a,kmin,'$s_t-q_t$','center','top')
plottext(b,kmin,'$s_t$','center','top')
stmp = 0.5*(a+b); ktmp = 0.5*f(stmp);
plottext(stmp,ktmp,'$C_t$','center','middle')

% Pot abandonment point
sstar = 0.25;
kstar = f(sstar);
plottext(sstar,kmin,'$s^*$','center','top')
plotvdash(sstar,kstar)
plothdash(sstar,kstar)
plottext(smin,kstar,'$p(0)$','right','middle')


%% FIGURE 7 Change in Consumer Surplus

% Initiate figure
figure
hold on

% Set plotting parameters
n = 1001;
qmin = 0;
qmax = 1;
pmin = 0;
pmax = 1;
p1 = 0.7;
p2 = 0.3;

% Define inverse demand curve
alpha = 0.15;
beta  = 1.25;
f = @(p) alpha*p.^(-beta);
q1 = f(p1);
q2 = f(p2);

% Plot area under inverse demand curve
p = nodeunif(n,0,pmax);
q = f(p);
for i=1:n
  if q(i)>q1 && p(i)>p2
    plot([0.005 f(p(i))],[p(i) p(i)],'y-','linewidth',2)
  end
end

% Plot inverse demand curve
plot(q,p)
xlim([qmin qmax])
ylim([pmin pmax])
set(gca,'xtick',[])
set(gca,'ytick',[])
axis square
plottext(qmin,pmin,'$0$','right','top')
plottext(qmax,pmin,'$q$','center','top')
plottext(qmin-0.02,pmax,'$p$','right','top')

% Annotate figure
plottext(qmin-0.02,p1,'$p_1$','right','middle')
plottext(qmin-0.02,p2,'$p_2$','right','middle')
plottext(q1,pmin,'$q_1$','center','top')
plottext(q2,pmin,'$q_2$','center','top')
plothdash(q1,p1)
plothdash(q2,p2)
plotvdash(q1,p1)
plotvdash(q2,p2)
plottext(0.9,0.36,'$p(q)$','right','top')

% To compute the change in consumer surplus `numerically'
[x,w] = qnwlege(15,p2,p1);
intn = w'*f(x);

% To compute the change in consumer surplus `analytically'
F = @(p) (alpha/(1-beta))*p.^(1-beta);
inta = F(p1)-F(p2);


%% SAVE FIGURES
printfigures(mfilename)