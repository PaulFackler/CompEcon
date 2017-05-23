% DEMFIN08 Affine asset pricing demonstration

disp(' ')
disp('DEMFIN08 Affine asset pricing demonstration')

% Define parameters
T     = 30;         % time-to-maturity
kappa = 0.1;        % speed of mean reversion
alpha = 0.05;       % long-run mean interest rate 
sigma = 0.1;        % interest rate volatility

% Convert parameters to standard form
a=kappa*alpha;
A=-kappa;
C=sigma;
b=0;
B=1;

% Call affine asset pricing solver
tau=linspace(0,T,301)';
[beta,beta0]=affasset(tau,a,A,b,B,C,1,0,0,0);

% Create plots
r=linspace(0,0.25,101)';
V=exp(beta0(end)+beta(end)*r);

figure(1)
plot(r,V); 
title([num2str(T) ' Year Zero-Coupon Bond Price'])
xlabel('r')
ylabel('V(r)')

figure(2)
plot(tau,[beta beta0]); 
title('\beta and \beta_0')
xlabel('Time to Maturity')

r=(0.03:0.01:0.08);
m=length(r);

warning off
R=-(beta0(:,ones(m,1))+beta*r)./tau(:,ones(m,1));
warning on
R(1,:)=r;
figure(3)
plot(tau,R)
xlabel('Time to Maturity')
ylabel('Yield')
title('Term Structures for Alternative Short Rates')

nn=length(r);
legend([repmat('r = ',nn,1) reshape(sprintf('%4.2f',r),4,nn)'])

prtfigs(mfilename,'Term Structures for Alternative Short Rates',[3])