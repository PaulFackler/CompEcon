% DEMBVP01 BVP illustrative example
function dembvp01
 disp(' ')
 disp('DEMBVP01 BVP illustrative example')

% Define the model
 A=[-1 -.5;0 -.5];
 model.func='pbvp01';
 model.tb=[0;1];
 model.params={A};

% Define the approximant
 n=5;
 a=0;
 b=2;
 fspace=fundefn('cheb',n-1,a,b);
 tnodes=funnode(fspace);
 fspace=fundefn('cheb',n,a,b);

% Initial conditions
 c=zeros(fspace.n,2);

% Evaluation points for plots
 tvals=linspace(a,b,201)';

% Call solver
 optset('broyden','defaults')
 [c,x,r]=bvpsolve(model,fspace,tnodes,c,tvals);

% Produce plots
 close all
 figure(1)
 C=exp(0.5)*(1-exp(-1));
 plot(tvals,exp(-tvals)-x(:,1),'-', ...
     tvals,C*exp(-tvals/2)+exp(-tvals)-x(:,2),'--');
 title('BVP Example: Approximation Errors')
 xlabel('t')
 ylabel('x-\phi(t)c')
 legend('x_1','x_2')

 figure(2)
 plot(tvals,r(:,1),'-',tvals,r(:,2),'--');
 title('BVP Example: Residual Functions')
 xlabel('t')
 ylabel('r')
 legend('x_1','x_2')
