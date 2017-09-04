%% Contents
%
% Demonstration programs to accompany
%   Lecture Notes in Computational Economic Dynamics
%   Professor Mario J. Miranda  
%   The Ohio State University
%   
% Based on book
%   Applied Computational Economics and Finance
%   Mario J. Miranda & Paul L. Fackler
%   2002, MIT Press, Cambridge MA
%
% These programs are continually reviseddem
%   This version 8/5/2016
%
%  Copyright(c) 1997-2016
%   Mario J. Miranda - miranda.4@osu.edu
%   Paul L. Fackler  - paul_fackler@ncsu.edu

% Introduction 
demintro01 % Inverse Demand Problem
demintro02 % Rational Expectations Agricultural Market Model
 
% Mathematics Review 
demmath01  % Taylor Approximations
demmath02  % Function Inner Products, Norms & Metrics
demmath03  % Discrete and Continuous Distributions
demmath04  % Standard Copulas
demmath05  % Inverse Function, Implicit Function and Mean Value Theorems
demmath06  % Operations with Markov Chains
 
% Linear Equations 
demlin01   % Linear Equations 
demlin02   % Ill-Conditioning of Vandermonde Matrices
demlin03   % Sparse Linear Equations
 
% Nonlinear Equations 
demslv01   % Root of f(x)=exp(-x)-1
demslv02   % Root of Rosencrantz Function
demslv03   % Fixedpoint of f(x) = x^0.5
demslv04   % Fixedpoint of g(x1,x2)= [x1^2+x2^3;x1*x2-0.5]
demslv05   % Newton and Broyden Paths
demslv06   % Nonlinear Equation Methods
demslv07   % Linear Complementarity Problem Methods
demslv08   % Nonlinear Complementarity Problem
demslv09   % Hard Nonlinear Complementarity Problem with Billup's Function
demslv10   % Linear Complementarity Problem
demslv11   % Rootfinding Reformulations of Nonlinear Complementarity Problem
demslv12   % Convergence Rates for Nonlinear Equation Methods
demslv13   % Spatial Equilibrium Model
 
% Finite-Dimensional Optimization
demopt01   % Maximization via Golden Search
demopt02   % Changes in Nelder-Mead Simplex
demopt03   % Nelder-Mead Simplex Method
demopt04   % Maximization of Rosencrantz Function by Various Methods
demopt05   % Optimization with qnewton
demopt06   % KKT Conditions for Constrained Optimization Problems
demopt07   % Bound-Constrained Optimization via Sequential LCP
demopt08   % Constrained Optimization with nlpsolve
 
% Quadrature
demqua01   % Equidistributed Sequences on Unit Square
demqua02   % Expectation of Function of Random Normal Vector
demqua03   % Area Under a Curve, Various Methods
demqua04   % Area under Normal PDF Using Simpson's Rule
demqua05   % Expected Utility Insurance Model
demqua06   % Area Under a Curve
demqua07   % Monte Carlo Simulation of Time Series

% Numerical Differentiation
demdif01   % Finite Difference Hessian Evaluation Structure
demdif02   % Error in Finite Difference Derivatives
demdif03   % Demonstrates fdjac and checkjac
demdif04   % Demonstrates fdhess
demdif05   % Finite-Difference Jacobians and Hessians
 
% Function Approximation
demapp01   % Approximating Functions on R
demapp02   % Approximating Functions on R^2
demapp03   % Standard Basis Functions and Nodes
demapp04   % Uniform- and Chebychev-Node Polynomial Approximation of Runge's Function
demapp05   % Chebychev Polynomial and Spline Approximation of Various Functions
demapp06   % Cournot Oligopoly Model
demapp07   % Compute Implicit Function via Collocation
demapp08   % Linear Spline Approximation
demapp09   % Monopolist's Effective Supply
 
% Discrete Time Discrete State Dynamic Programming
demddp01   % Mine Management Model
demddp02   % Asset Replacement Model
demddp03   % Asset Replacement Model With Maintenance
demddp04   % Binomial American Put Option Model
demddp05   % Water Management Model
demddp06   % Bioeconomic Model
demddp07   % Renewable Resource Model
demddp08   % Job Search Model
demddp09   % Deterministic Cow Replacement Model
demddp10   % Stochastic Cow Replacement Model
demddp11   % Stochastic Optimal Growth Model

% Discrete Time Continuous State Dynamic Programming
demdp01a   % Timber Harvesting - Direct Solution, Linear Approximation
demdp01b   % Timber Harvesting - Direct Solution, Cubic Spline Approximation
demdp01c   % Timber Harvesting - DPSOLVE Solution, Cubic Spline Approximation
demdp02    % Asset Replacement Model
demdp03    % Industry Entry-Exit Model
demdp04    % Job Search Model
demdp05    % American Put Option Pricing Model
demdp06    % Ramsey Economic Growth Model
demdp07    % Stochastic Economic Growth Model
demdp08    % Public Renewable Management Model
demdp09    % Private Non-Renewable Resource Management Model
demdp10    % Water Resource Management Model
demdp11    % Monetary Policy Model
demdp12    % Production Management Model
demdp13    % Inventory Management Model
demdp14    % Livestock Feeding Model (Euler Conditions)
demdp15    % Saving with Transactions Costs Model
demdp16    % Linear-Quadratic Model
demdp17    % Miscellaneous Lecture Note Figures
demdp18    % Credit With Technology Adoption Model
demdp19    % Credit with Strategic Default Model (Broyden)
demdp20    % Lifecycle Consumption-Savings Model
demdp21    % Savings and Insurance Model - Base Case
demdp22    % Savings and Insurance Model - Sensitivity Analysis

% Rational Expectations Models
demrem01   % Asset Pricing Model
demrem02   % Commodity Storage Model
demrem03   % Government Price Support Model
 
% Dynamic Games
demgame01  % Production Capacity Game Model
demgame02  % Income Redistribution Game Model
demgame03  % Marketing Board Game Model

% Ordinary Differential Equations
demode01   % Stability of Linear ODEs
demode02   % Generic Nonlinear ODE Example
demode03   % IVP Linear ODE Example
demode04   % Non-IVP Linear ODE Example
demode05   % Commodity Storage Model
demode06   % Predator-Prey Model
demode07   % Commercial Fisheries Model
demode08   % Lorentz Strange Attractor
demode09   % Tobin's Q
demode10   % Regional Migration Model

% Continuous Time Deterministic Optimal Control
demdoc01   % Deterministic Optimal Consumption-Investment Model
demdoc02   % Deterministic Optimal Economic Growth Model
demdoc03   % Deterministic Nonrenewable Resource Model
demdoc04   % Deterministic Renewable Resource Model
demdoc05   % Deterministic Production Adjustment Model

% Continuous Time Stochastic Optimal Control
demsoc00   % Ito Processes
demsoc01   % Stochastic Consumption-Investment Model
demsoc02   % Stochastic Portfolio Model
demsoc03   % Stochastic Optimal Economic Growth Model
demsoc04   % Stochastic Renewable Resource Model
demsoc05   % Stochastic Optimal Fish Harvest Model
demsoc06   % Stochastic Production-Adjustment Model
demsoc07   % Stochastic Nonrenewable Resource Model

% Regime Switching Models
demrs01    % Asset Abandonment Model
demrs02    % Fish Harvest Model
demrs03    % Dixit Entry/Exit Model
 
% Impulse Control Models
demic01    % Asset Replacement Model
demic02    % Timber Harvesting Model
demic03    % Storage Management Model
demic04    % Capacity Choice Model
demic05    % Cash Management Model
demic06    % Fish Harvest Model