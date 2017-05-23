function [out1,out2,out3] = pgame03(flag,s,x,e,smax,eps,gamma,A,beta,theta,alpha,cost,xi);
% Function file for capital investment game

n = size(s,1);
ds = 2;
dx = 2;

switch flag
case 'b'; % BOUND FUNCTION
  xl = zeros(n,dx);
  xu = zeros(n,dx);
  xu = smax(ones(n,1),:) - (1-xi)*s;
  out1=xl; out2=xu; out3=[];
case 'f1'; % REWARD FUNCTION
%  f   = zeros(n,1);
%  fx  = zeros(n,dx);
%  fxx = zeros(n,dx,dx);
  Prof = Profit(s,n,alpha,beta,gamma,A,eps);
  c = cost(:,1);
  xx = x(:,1);
  if theta == 1;
    f =  Prof(:,1) - c*xx;
    fx =  - c*ones(n,1);
    fxx = zeros(n,1);
  else
    f = Prof(:,1) - c*xx.^theta/theta;
    fx =  - c*xx.^(theta-1);
    fxx = - (theta - 1)*c*xx.^(theta - 2);
  end;
  out1 = zeros(n,1);
  out2 = zeros(n,dx);
  out3 = zeros(n,dx,dx);
  out1 = f;
  out2(:,1)= fx;
  out3(:,1,1) = fxx;
% out1=f; out2=fx; out3=fxx;
case 'f2'; % REWARD FUNCTION
%  f   = zeros(n,1);
%  fx  = zeros(n,dx);
%  fxx = zeros(n,dx,dx);
  Prof = Profit(s,n,alpha,beta,gamma,A,eps);
  c = cost(:,2);
  xx = x(:,2);
  if theta == 1;
    f =  Prof(:,2) - c*xx;
    fx =  - c*ones(n,1);
    fxx = zeros(n,1);
  else
    f = Prof(:,2) - c*xx.^theta/theta;
    fx =  - c*xx.^(theta-1);
    fxx = - (theta - 1)*c*xx.^(theta - 2);
  end;
  out1 = zeros(n,1);
  out2 = zeros(n,dx);
  out3 = zeros(n,dx,dx);
  out1 = f;
  out2(:,2)= fx;
  out3(:,2,2) = fxx;
% out1=f; out2=fx; out3=fxx;
case 'g'; % STATE TRANSITION FUNCTION
  g   = zeros(n,ds);
  gx  = zeros(n,ds,dx);
  gxx = zeros(n,ds,dx,dx);
  g = (1-xi)*s + x;
  gx(:,1,1) = ones(n,1);
  gx(:,2,2) = ones(n,1);
  out1=g; out2=gx; out3=gxx;
otherwise; error('Improper function flag')
end 


function Prof = Profit(s,n,alpha,beta,gamma,A,eps);

  Q = ones(2,n);
  aalpha = alpha'*ones(1,n);
  bbeta = beta'*ones(1,n);
  ggamma = gamma'*ones(1,n);
  AA = A'*ones(1,n);
  ee = (1 + diag(eps))*ones(1,n);
  eeps = eps + diag(1 - alpha);
  RHS = log(ggamma.*aalpha) + bbeta.*log(s') - log(AA.*ee);
  Q = exp(eeps\RHS);
  P = exp(log(AA) + eps*log(Q));
  C = exp(log(ggamma) + aalpha.*log(Q) + bbeta.*log(s'));
  Prof = (P.*Q - C)';

