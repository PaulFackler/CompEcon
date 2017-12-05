% COMPECON: COMPUTATIONAL ECONOMICS AND FINANCE
% Toolbox functions to accompany:
%   Applied Computational Economics and Finance
%   Mario J. Miranda & Paul L. Fackler
%   2002, MIT Press, Cambridge MA
%   ISBN 0-262-13420-9
%
% ROOTFINDING AND OPTIMIZATION
%   BISECT    Uses method of bisection to find roots for 1-D functions
%   BROYDEN   Computes root of function via Broyden's Inverse Method
%   BROYDENX  Computes root of function via Broyden's Method
%   FIXPOINT  Computes fixpoint of function using function iteration
%   GJACOBI   Solves Ax=b using Jacobi iteration
%   GOLDEN    Computes local maximum of univariate function on interval via Golden Search
%   GOLDENX   Computes local maximum of univariate function on interval via Golden Search
%   GSEIDEL   Solves Ax=b using Gauss-Seidel iteration
%   LCPBAARD  Uses block Baard method to solve linear complementarity problem
%   LCPLEMKE  Solves linear complementarity problem using Lemke's algorithm
%   LCPSOLVE  Solves linear complementarity problem using safeguarded Newton method
%   LEMKE     Solves linear complementarity problems (LCPs).
%   LPSOLVE   Solves linear programming problems
%   MINMAX    Minimax transformation for solving NCP as rootfinding problem
%   NCPJOSE   Solves nonlinear complementarity problem using sequential LCP method
%   NCPSOLVE  Solves nonlinear complementarity problem
%   NELDMEAD  Maximizes function via Nelder-Mead algorithm
%   NEWTON    Computes root of function via Newton's Method with backstepping
%   OPTSTEP   Solves a one dimensional optimal step length problem
%   QNEWTON   Solves unconstrained maximization problem using quasi-Newton
%   SMOOTH    Reformulates an MCP as a semismooth function
%
% QUADRATURE
%   QNWBETA   Computes quadrature nodes and weights for Beta(a,b) distribution
%   QNWCHEB   Computes multivariate Guass-Chebyshev quadrature nodes and weights
%   QNWEQUI   Generates equidistributed sequences
%   QNWGAMMA  Quadrature nodes and weights for Gamma(a) distribution
%   QNWLEGE   Computes multivariate Guass-Legendre quadrature nodes and weights
%   QNWLOGN   Computes Gauss-Hermite nodes and weights multivariate lognormal distribution
%   QNWNORM   Computes nodes and weights for multivariate normal distribution
%   QNWSIMP   Computes multivariate Simpson quadrature nodes and weights
%   QNWTRAP   Computes multivariate trapezoid rule quadrature nodes and weights
%   QNWUNIF   Computes nodes and weights for multivariate uniform distribution
%   QUADRECT  Integrates function on a rectangular region in R^n
%
% FUNCTION APPROXIMATION
%   CHEBBAS   Computes basis matrices for Chebyshev polynomials
%   CHEBDEF   Defines parameters for Chebyshev polynomial functions
%   CHEBDOP   Creates differential operator matrices for Chebyshev polynomials.
%   CHEBNODE  Computes standard nodes for Chebyshev polynomials
%   FOURBAS   Defines basis matrices for Fourier series
%   FOURDEF   Defines parameters for Fourier basis functions
%   FOURDOP   Computes differential operator for Fourier functions
%   FOURNODE  Computes standard nodes for Fourier basis
%   FUNBAS    Computes a basis matrix
%   FUNBASX   Creates basis structures for function evaluation
%   FUNBCONV  Converts among basis structure formats
%   FUNCONV   Converts from one basis family to another
%   FUND      Evaluates functions and first 2 derivatives
%   FUNDEF    Creates a fauction family definition structure and/or performs checks
%   FUNDEFN   Defines a function family structure
%   FUNDOP    Computes derivative operators
%   FUNEVAL   Evaluates multivariate functions with linear bases.
%   FUNFITF   Computes interpolation coefficients for D-dim function.
%   FUNFITXY  Computes interpolation coefficients for d-dim function.
%   FUNHESS   Computes the Hessian for FUN functions 
%   FUNJAC    Computes the Jacobian for FUN functions 
%   FUNNODE   Computes default nodes for a family of functions
%   LINBAS    Piecewise linear basis functions
%   LINDEF    Computes standard breakpoints for linear spline
%   LINDOP    Differential operator for a piecewise linear function
%   LINNODE   Standard nodes for linear spline
%   SPLIBAS   Computes polynomial spline basis.
%   SPLIDEF   Defines default parameters for spline functions
%   SPLIDOP   Creates differential operator matrices for polynomial splines.
%   SPLINODE  Computes standard nodes for splines using knot averaging.
%
% DISCRETE TIME DYNAMIC MODELING
%   DDPSIMUL   Monte Carlo simulation of discrete-state/action controlled Markov process
%   DDPSOLVE   Solves discrete-state/action dynamic program
%   DPCHECK    Checks derivatives for dp files
%   DISCRAND   Discrete random variable simulator
%   DPSIMUL    Monte Carlo simulation of discrete time controlled Markov process
%   DPSOLVE    Solves discrete time continuous-state/action dynamic program
%   DPSTST     Computes invariant distribution for continuous-state/action controlled dynamic program
%   GAMESOLVE  Solves discrete time continuous-state/action Bellman equations for dynamic games
%   LQAPPROX   Forms and solves linear-quadratic approximation of DP model
%   MARKOV     Analyzes Markov transition probability matrices
%   REMSIMUL   Simulates state paths in rational expectations models
%   REMSOLVE   Solves rational expectations models
%   REMSTST    Computes invariant distribution for rational expectations models
% 
% ODE SOLVERS AND CONTINUOUS TIME DYNAMIC MODELLING
%   AFFASSET      Solves affine asset pricing models
%   BVPSOLVE      Solves general first order boundary value problems
%   CTBASEMAKE    Basis matrices for continuous time collocation
%   CTSTEADYSTATE Finds deterministic steady state for continuous time models
%   FINDSTATE     Calibrates an asset pricing model to data
%   FINSOLVE      Solves continuous time asset pricing problems
%   ICSOLVE       Solves continuous time impulse control models
%   ITODENSITY    Long-run densities for 1-D Ito processes
%   ITOSIMUL      Monte Carlo simulation of a (possibly) controlled Ito process
%   RK4           Solves initial value problems using fourth-order Runge-Kutta
%   RSSOLVE       Solves continuous time regime switching models
%   SCSOLVE       Solves stochastic control problems
%
% UTILITIES
%   CHECKJAC  Compares analytic and finite difference derivative
%   CHKFIELDS Checks if a variable S is a valid structure with fields F
%   CKRON     Repeated Kronecker products on a cell array of matrices
%   CKRONX    The product of repeated Kronecker products and a matrix
%   CKRONXI   The product of repeated inverse Kronecker products and a matrix
%   CSIZE     Returns dimension information for cell arrays
%   GETINDEX  Finds the index value of a point
%   GRIDMAKE  Forms grid points
%   INDEX     Converts between single and multiple indices
%   KERNEL    Computes a kernel estimate of a PDF
%   TABLOOKUP Performs a table lookup
%   MEXALL    Creates MEX files for CompEcon toolbox
%   MINTERP   Multidimensional interpoation
%   NODEUNIF  Computes uniform nodes for intervals R^n
%   OPTGET    Utility to get previously set function default values 
%   OPTSET    Utility to set function options 
%
% SPECIAL FUNCTIONS AND MISC.
%   BAW       Barone-Adesi/Whaley American option pricing model
%   BS        Black-Scholes option pricing model
%   CDFN      Computes the CDF of the standard normal distribution
%   CONFHYP   Computes the confluent hypergeometric function 
%   DIGAMMA   Computes the digamma (psi) function for positive arguments
%   FDHESS    Computes finite difference Hessian
%   FDJAC     Computes two-sided finite difference Jacobian
%   FHESS     Alternative finite difference Hessian procedure
%   FJAC      Alternative finite difference Jacobian procedure
%   IMPVOL    Computes option implied volatilities
%   MONTNORM  Computes pseudo-random multivariate normal variates
%   PSI       Calculates the value of the psi (digamma) function
%   TRIGAMMA  Calculates the value of the trigamma function