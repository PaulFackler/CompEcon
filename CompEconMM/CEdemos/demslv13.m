%% DEMSLV13 Spatial Equilibrum Model
%
% Solves three-country spatial equilibrium with linear demands and supplies
% as linear complementarity problem.

function demslv13

% Preliminary tasks
demosetup(mfilename)

% Model parameters
as = [9; 3; 18];
bs = [1; 2; 1];
ad = [42; 54; 51];
bd = [3; 2; 1];
c = [0 3 9;3 0 3;6 3 0];
a = zeros(9,1);
b = inf*ones(9,1);

% Solve spatial equilibrium model
optset('ncpsolve','ShowIters',1)
x = zeros(9,1);
x = ncpsolve(@f,a,b,x,as,bs,ad,bd,c);

% Output
x = reshape(x,3,3);
p = as + bs.*sum(x,2);

fprintf('\n\nEquilibrium Trade Flows\n')
fprintf('From\\To   Country 1   Country 2   Country 3\n')
fprintf('Country 1   %5.1f       %5.1f       %5.1f\n',x(1,:))
fprintf('Country 2   %5.1f       %5.1f       %5.1f\n',x(2,:))
fprintf('Country 3   %5.1f       %5.1f       %5.1f\n',x(3,:))

fprintf('\n\nEquilibrium Prices\n')
fprintf('Country 1   %5.1f\n',p(1))
fprintf('Country 2   %5.1f\n',p(2))
fprintf('Country 3   %5.1f\n',p(3))


function [fval,fjac] = f(x,as,bs,ad,bd,c)
x = reshape(x,3,3);
ps = as + bs.*(sum(x,2));
pd = ad - bd.*(sum(x,1))';
ps = ps(:,ones(1,3));
pd = pd(:,ones(1,3))';
fval = pd - ps - c;
fval = reshape(fval,9,1);
fjac = [];