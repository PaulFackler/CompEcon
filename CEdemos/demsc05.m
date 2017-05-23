% DEMSC05 Sequential learning Model
 disp(' ')
 disp('DEMSC05 Sequential learning Model')

% ************* Initialize Parameters ****************
  r=0.05; 
  delta=0.05; 
  sigma=.2;
  Qm=20;
  cbar=10;
  C=40;
  N=300;
  n=15;
  hi=log(100);

% ********* Call solution function ********* 
  [Q,Pstar,c,cdef,A1,A2,beta1,beta2] = ...
       learn(r,delta,sigma,Qm,cbar,C,N,n,hi);

% ********* Produce plots **********

% Compute deterministic free boundary
  gam=log(C/cbar)/Qm;
  pstar0=(r*exp(gam*(Qm-Q))+gam*exp(-r*(Qm-Q)))/(r+gam)*cbar;

  close all
  figure(1)
  plot(Q,Pstar,'k',[Qm;Qm],[0;5*cbar],'k:',[Qm;1.4*Qm],[cbar;cbar],'k')
 % Plot bounds on free boundary
%  hold on; plot(Q,[pstar0;cbar*exp(gam*(Qm-Q))],'k--'); hold off 
  title('Optimal Activation Boundary');
  xlabel('Q')
  ylabel('P')
  axis([0 1.4*Qm 0 5*cbar])
  TextFontSize=8;
  h=text(5,40,'Production Region:');
  set(h,'FontSize',TextFontSize);
  h=text(6,37,'V(P,Q) must be');
  set(h,'FontSize',TextFontSize);
  h=text(7,34,'computed numerically');
  set(h,'FontSize',TextFontSize);
  h=text(3,11,'Non-Production Region:');
  set(h,'FontSize',TextFontSize);
  h=text(4,8,'A(Q) computed from');
  set(h,'FontSize',TextFontSize);
  h=text(5,5,'value matching condition');
  set(h,'FontSize',TextFontSize);
  h=text(20.5,38,'Production');
  set(h,'FontSize',TextFontSize);
  h=text(23,35,'Region:');
  set(h,'FontSize',TextFontSize);
  h=text(22,32,'V known');
  set(h,'FontSize',TextFontSize);
  h=text(20.5,8,'Non-Production');
  set(h,'FontSize',TextFontSize);
  h=text(23,5,'Region:');
  set(h,'FontSize',TextFontSize);
  h=text(22,2,'V known');
  set(h,'FontSize',TextFontSize);
  h=text(20.5,48,'Q_m');
  set(h,'FontSize',TextFontSize);
  h=text(3,25,'P^*');
  set(h,'FontSize',TextFontSize);

% Plot value functions V(P,Q_k) for several Q_k.
  figure(2)
  Y=linspace(0,hi,101)';
  ind=lookup(Q,[0;5;10;15;20]);  % interpolate at several values of Q
  p1=exp(Y)*Pstar(ind)';
  V1=funeval(c(1:n,ind),cdef,Y);
  v0=funeval(c,cdef,0);
  A1=v0.*(Pstar.^(-beta1))';
  p0=linspace(0,1,101)'*Pstar';
  V0=A1(ones(101,1),:).*(p0.^beta1);
  plot(p1,V1,'k',p0(:,ind),V0(:,ind),'k',Pstar(ind),V1(1,:),'k*')
  title('Value Function');
  xlabel('P')
  ylabel('V(P,Q)')
  xlim([0,4*cbar])
  ylim([0,500])
  text(12,300,'Q=20=Q_m')
  text(30,200,'Q=0')

disp('Critical price for stochastic and deterministic cases (Q=20)')
disp([Pstar(1) pstar0(1)])

prtfigs(mfilename,'Solution to the Sequential Learning Model',[1 2])