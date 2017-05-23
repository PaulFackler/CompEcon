% Demonstration Files for CompEcon Toolbox
% * indicates a demo that produces figures in the text

% Default settings for printing
%set(0,'DefaultAxesPosition',[0.1    0.25    0.55    0.6667]);
%set(0,'DefaultFigurePaperPosition',[1.5 3.625 5.5 3.85]);

set(0,'DefaultFigurePaperPosition',[4.25 3.3125 2.5 1.875]);

set(0,'DefaultAxesFontSize',8);
set(0,'DefaultAxesFontName','Helvetica');


optset('prtfigs','dir','d:\booknew\figs\')
optset('savefigs','dir','d:\booknew\savefigs\')

set(0,'DefaultAxesPosition',[0.15  0.15  0.7  0.75]);
set(0,'DefaultFigurePaperPosition',[1.5 3.625 5.5 3.85]);
set(0,'DefaultFigureClipping','on');
set(0,'DefaultAxesFontName','Times-Roman');
set(0,'DefaultTextFontName','Times-Roman')
set(0,'DefaultAxesFontSize',16);
set(0,'DefaultTextFontSize',16);

% Rootfinding Solvers
demslv11 % Compares fixed-point, Newton & Secant iterations *
demslv12 % Cournot oligopoly example *
demslv13 % Illustrates CP problem *
demslv14 % Illustrates minmax reformulation of CP problems *
demslv15 % Illustrates minmax and semi-smooth reformulations of CP problems *
demslv16 % A Difficult NCP *

% 
% Optimization Solvers
demopt01 % Demonstrates golden search method *
demopt02 % Displays changes in a simplex *
demopt03 % Demonstrates Nelder-Mead simplex method *
demopt04 % Maximization of banana function via Quasi-Newton *
demopt05 % Illustrates step length algorithms *
demopt06 % Illustrates relationship between maximization and CPs *

% 
% Quadrature Methods 
demqua01 % Alternative Equidistributed Sequences *
% 
% Derivatives and Differential Equation Examples
demdif01 % Creates a plot to illustrate finite difference evaluation *
demdif02 % Errors in Numerical Derivatives *
demdif03 % Commercial fisheries example *
% 
% Function Approximation
demapp01  % Plots Basis Functions and standard nodes *
demapp02  % Compares Alternative Approximates for exp(-x) *
demapp04  % Compares Alternative Approximates for Runge's function *
demapp05  % Compares Alternative Approximates for several functions *
demapp06  % Approximate exp(-x) using polynomials and cubic splines *
demapp09  % Cournot oligopoly example *
demapp11  % Equilibrium Storage Problem
% 
% Discrete Dynamic Programming
demddp01   % Mine Management Problem
demddp02   % Machine Replacement Model
demddp03   % Asset Replacement-Maintenance Model
demddp04   % Binomial Americal Option Pricing Model
demddp06   % Water Management Model
demddp07   % Bioeconomic Model
