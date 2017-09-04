% RUNDEMOS Runs all demonstration files the accompany the CompEcon Toolbox
% 
% COMPECON: COMPUTATIONAL ECONOMICS AND FINANCE
% Demonstration files to accompany:
%   Applied Computational Economics and Finance
%   Mario J. Miranda & Paul L. Fackler
%   2002, MIT Press, Cambridge MA
%   ISBN 0-262-13420-9


% *********** Introductory problems *********** 
demint01 % Solves inverse demand problem *
demint02 % Solves nonlinear rational expectations agricultural market model
%
% *********** Linear Equations *********** 
demlin01 % Compares effort to solve Ax=b by different methods
demlin02 % Demonstrates ill-conditionsing of Vandermonde matrix
demlin03 % Compares various methods for solving sparse linear equations
%
% *********** Rootfinding Problems *********** 
demslv01 % Compute root of Rosencrantz function via Newton and Broyden methods
demslv02 % Computes root of function using Newton, quasi-Newton and Bisection
demslv03 % Computes fixed point y=[x1^2+x2^3;x1*x2-0.5] via function iteration and Newton methods
demslv04 % Graphical demonstration of computing fixed point y=x^0.5
demslv05 % Compare various linear complementarity algorithms
demslv06 % Demonstrates different NCP methods
demslv07 % Demonstrates different NCP methods with Billup's function
demslv08 % Compute fixedpoint of using function iteration and Newton method
demslv10 % Demonstrates bisection method
demslv11 % Generates figures for Chapter 3
demslv12 % Cournot demonstration
demslv13 % Illustrates CP problems
demslv14 % Illustrates Problematic Complementarity Problems
demslv15 % Illustrates minmax and semi-smooth reformulations of CP problems
demslv16 % Billups Hard NCP
% 
% *********** Optimization Problems ***********
demopt01 % Illustrates maximization via golden search
demopt02 % Displays changes in a simplex
demopt03 % Demonstrates Nelder-Mead simplex method
demopt04 % Demonstrates Quasi-Newton maximization
demopt05 % Step length determination
demopt06 % Illustrates constrained optimization problems
demopt07 % Maximizes function subject to simple bounds via sequential LCP
% 
% *********** Quadrature Methods ***********
demqua01 % Plots equi-distributed sequences in 2-D
demqua02 % Compute expectation of function of random normal vector
demqua03 % Compares quadrature methods
demqua04 % Compare various quadrature methods
demqua05 % Compares Chebyshev and Legendre quadrature nodes and weights
% 
% *********** Derivatives and Differential Equation Examples ***********
demdif01 % Plot to illustrate finite difference Hessian evaluation
demdif02 % Demonstrates errors in finite difference differentiation of e(x)
demdif03 % Commercial Fishery Model (from V.L. Smith)
% 
% *********** Function Approximation ***********
demapp00 % Creates Table 6.1
demapp01 % Plot basis functions and standard nodes for major approximation schemes
demapp02 % Compare polynomial and spline approximation of function exp(-x)
demapp03 % Compare conditioning of Vandermonde and Chebychev matrices
demapp04 % Approximating Runge's function and other comparisons
demapp05 % Demonstrates alternative approximants for various functions
demapp06 % Approximate exp(-x) on [-1,1] via spline and Chebychev approximation
demapp07 % Approximate y=x1/exp(x2) on [0,5]x[-1,1] using spline and Chebychev approximation
demapp08 % Approximate y=(exp(x1*x2)+x3^3)/100 on [0,1]^3 via spline and Chebychev approximation
demapp09 % Cournot oligopolist problem
demapp10 % Approximate function inverse via collocation
demapp11 % Equilibrium Storage Problem
% 
% *********** Boundary Value Problems ***********
dembvp01 % BVP illustrative example 
dembvp02 % Commodity storage example

% 
% *********** Discrete Dynamic Programming ***********
demddp01 % Mine Management Model
demddp02 % Asset Replacement Model
demddp03 % Asset Replacement with Maintenance Model
demddp04 % Binomial Americal Option Pricing Model
demddp05 % Water Management Model
demddp06 % Bioeconomic Model
demddp07 % Renewable Resource Model 
demddp08 % Job Search Model
demddp09 % Deterministic Cow Replacement Model
demddp10 % Stochastic Cow Replacement Model
demddp11 % Optimal growth model
%
% *********** Dynamic Programming ***********
demdp01 % Asset Replacement Model
demdp02 % Industry Entry-Exit Model
demdp03 % Timber Cutting Model
demdp04 % American Option Pricing Model
demdp05 % Job Search Model
demdp06 % Asset Replacement-Maintanence Model
demdp07 % Optimal Growth Model
demdp08 % Renewable Resource Model
demdp09 % Non-Renewable Resource Model
demdp10 % Water Management Model
demdp11 % Monetary Policy Model
demdp12 % Production-Adjustment Model
demdp13 % Producton-Inventory Model
demdp14 % Livestock Feeding Model
%
% *********** Dynamic Games ***********
if 0 %  Change to run the games models (these take a long time)
demgame01 % Capital Production Game
demgame02 % Income Redistribution Game
demgame03 % Marketing Board Game
end
% 
% *********** Rational Expectations Models ***********
demrem01 % Asset Pricing Model
demrem02 % Commodity Storage Model
demrem03 % Government Price Support Model
%  
% *********** Financial Asset Pricing ***********
demfin01 % Cox-Ingersoll-Ross bond pricing example
demfin02 % Black-Scholes option pricing model
demfin03 % Heston's stochastic volatility option pricing model
demfin04 % American put option pricing problem
demfin05 % Barrier option pricing problem
demfin06 % Compound option pricing example - option on a bond
demfin07 % Asian option pricing demonstration
demfin08 % Affine asset pricing demonstration
demfin09 % Financial asset calibration demonstration
% 
% *********** Stochastic Control ***********
demsc01 % Optimal Growth Model
demsc02 % Renewable Resource Model
demsc03 % Production-Adjustment Model
demsc04 % Optimal Fish Harvest Model
demsc05 % Sequential Learning Model
demsc06 % Nonrenewable Resource Model (continuous time)
% 
% *********** Regime Switching Problems ***********
demrs01 % Asset Abandonment Model
demrs02 % Optimal Fish Harvest Model
demrs03 % Entry/exit Model
%
% *********** Impulse Control Problems ***********
demic01 % Asset Replacement Demonstration
demic02 % Timber Harvesting Demonstration
demic03 % Storage Management Demonstration
demic04 % Capacity Choice Demonstration
demic05 % Cash Management Demonstration
demic06 % Optimal Fish Harvest problem (infinite harvest rate)