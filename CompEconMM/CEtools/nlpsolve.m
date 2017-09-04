%% NLPSOLVE   
%
%  Solves nonlinear constrained maximization problem max f(x) s.t. x>=0, g(x)<=b
%
%  Think of
%   x as levels of n economic activities 
%   b as amounts available of m resources
%   f(x) as profit from the n activities
%   g(x) as amount of resources required to sustain activity levels x
%
%  Usage
%   [x,fval,lambda,MP,exitflag] = nlpsolve(x,f,g,b)
%  Let
%    n        : dimension of domain of f and g, number of activities
%    m        : dimension of range of g, number of resources
%  Input
%    x        : n.1 initial guess for constrained maximum
%    f        : function from R^n to R of form fval=f(x)
%    g        : function from R^n to R^m of form gval=g(x)
%    b        : m.1 vector
%  Output
%    x        : n.1 constrained maximum
%    fval     : optimal value of objective
%    lambda   : m.1 shadow prices of resources
%    MP       : n.1 economic marginal profits of activities
%    exitflag : 1 if otimality conditions satisfied, otherwise see fmincon
%  Note
%    This function simply reformats data for canonical constrained
%    maximization problems to conform to the usage of the Matlab
%    constrained minimization routine fmincon.

%  Copyright(c) 1997-2015
%   Mario J. Miranda - miranda.4@osu.edu
%   Paul L. Fackler  - paul_fackler@ncsu.edu

function [x,fval,lambda,MP,exitflag] = nlpsolve(x,ff,gg,bb)

clear f b g
global f g b
f = ff;
g = gg;
b = bb;
warning('off')
options = optimoptions('fmincon','Display','none','TolFun',1e-10,'TolX',1e-10,'FinDiffType','central');
[x,fval,exitflag,output,lambda,grad,hessian] = fmincon(@F,x,[],[],[],[],[0;0],[],@C,options);
fval - -fval;
MP = -lambda.lower;
lambda = lambda.ineqnonlin;
disp(['  '])
disp(['Maximum Objective'])
disp([f(x)])
disp(['         x        MP'])
disp([x MP])
disp(['    lambda   surplus'])
disp([lambda b-g(x)])
warning('on')

function y = F(x)
global f
y = -f(x);

function [c,ceq] = C(x)
global g b
ceq = [];
c = g(x)-b;