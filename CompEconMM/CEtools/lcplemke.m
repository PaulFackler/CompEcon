%% LCPLEMKE
%
%  Solves linear complementarity problem using Lemke's algorithm
%
%  Problem
%     z = M*x + q
%     a <= x <= b
%     x_i > a_i => z_i => 0
%     x_i < b_i => z_i =< 0
%
%  Usage
%    [x,z,err] = lcplemke(M,q,a,b,x)
%  Input
%    M         : n.n matrix
%    q         : n.1 vector
%    a         : n.1 lower bound on x
%    b         : n.1 upper bound on x
%    x         : n.1 initial guess for solution
%  Output
%    x         : solution to lcp
%    z         : function value at x
%    err       : 0-Solution found
%                1-Maximum iterations exceeded
%                2-Unbounded ray termination
%                3-Initial basis cannot be computed - try new value of x
%  Options
%    maxit     : maximum number of iterations (min(1000,25*n))
%    zer_tol   : convergence tolerance (1e-5)
%    piv_tol   : convergence tolerance (1e-8)

%  Copyright(c) 1997-2015
%   Mario J. Miranda - miranda.4@osu.edu
%   Paul L. Fackler  - paul_fackler@ncsu.edu

function [x,z,err] = lcplemke(M,q,a,b)

n = length(q);
if size(M)~=[n n]
  error('Matrices are not compatible');
end

% Set option defaults, if not set by user with OPTSET
maxit    = optget('lcplemke','maxit',min(1000,25*n));
zer_tol  = optget('lcplemke','zer_tol',1e-5);
piv_tol  = optget('lcplemke','piv_tol',1e-8);

Minv = inv(M);
MM = [-Minv Minv ; Minv -Minv];
z0 = [ max(-M*q-q,0) ; max(M*q+q,0) ];
qq = [-M\q-a ; M\q+b];
n  = 2*n;
err = 0;

% Trivial solution exists
if all(qq>=0)
  z = zeros(n,1); return;
end

% Initializations 
z = zeros(2*n,1);
iter = 0;
x = qq;
t = 2*n+1;

% Determine initial basis
if isempty(z0)
  bas = [];
  nonbas = (1:n)';
else
  bas = find(z0>0);
  nonbas = find(z0<=0);
end

B = -speye(n);

% Determine initial values
if ~isempty(bas)
  B = [MM(:,bas) B(:,nonbas)];
  if condest(B)>1e16
    z = []; 
    err = 3; 
    return
  end
  x = -(B\qq);
end

% Check if initial basis provides solution
if all(x>=0)
  z(bas) = x(1:length(bas));
  z = z(1:n);
  return
end

% Determine initial leaving variable
[tval,lvindex] = max(-x);
bas = [bas;(n+nonbas)];
leaving = bas(lvindex);

bas(lvindex) = t;                     % pivot in the artificial variable

U = x<0;
Be = -(B*U);
x = x+tval*U;
x(lvindex) = tval;
B(:,lvindex) = Be;

% Main iterations begin here
for iter=1:maxit
  % Check if done; if not, get new entering variable
  if (leaving == t)
    break
  elseif (leaving <= n)
    entering = n+leaving;
    Be = zeros(n,1); 
    Be(leaving) = -1;
  else
    entering = leaving - n;
    Be = MM(:,entering);
  end
  d = B\Be;
  
  % Find new leaving variable
  j = find(d>piv_tol);                % indices of d>0
  if isempty(j)                       % no new pivots - ray termination
    err = 2;
    break
  end
  theta = min((x(j)+zer_tol)./d(j));  % minimal ratios, d>0
  j = j((x(j)./d(j))<=theta);         % indices of minimal ratios, d>0
  lvindex = find(bas(j)==t);          % check if artificial among these
  if ~isempty(lvindex)                % always use artifical if possible
    lvindex = j(lvindex);
  else                                % otherwise pick among set of max d
    theta = max(d(j));
    lvindex = find(d(j)==theta);
    lvindex = ceil(length(lvindex)*rand(1,1));  % if multiple choose randomly
    lvindex = j(lvindex);
  end
  leaving = bas(lvindex);
  
  % Perform pivot
  ratio = x(lvindex)./d(lvindex);
  x = x - ratio*d;
  x(lvindex) = ratio;
  B(:,lvindex) = Be;
  bas(lvindex) = entering;
end                                   % end of iterations
if iter>=maxit && leaving~=t
  err = 1;
end

z(bas) = x; z = z(1:n);

% Display warning messages if no error code is returned
if nargout<2 && err~=0
  s = 'Warning: solution not found - ';
  if err==1
    disp([s 'Iterations exceeded limit']);
  elseif err==2
    disp([s 'Unbounded ray']);
  elseif err==3
    disp([s 'Initial basis infeasible']);
  end
end

z = z(n/2+1:n)-z(1:n/2);
x = M\(z-q);