% DEMAPP06 Approximate exp(-x) on [-1,1] via spline and Chebychev approximation
disp('DEMAPP06 Approximate exp(-x) on [-1,1] via spline and Chebychev approximation')
close all

% Set endpoints of interpolation interval
a =  -1;                            % left endpoint
b =   1;                            % right endpoint
n =  10;                            % order of interpolatioin

% Construct refined uniform grid for error ploting
x = nodeunif(1001,a,b);

% Compute actual and fitted values on grid
[y,d,s] = fapp06(x);                % actual

% Construct and evaluate Chebychev interpolant
chebproj = fundefn('cheb',n,a,b);   % construct projection space
c     = funfitf(chebproj,'fapp06'); % project onto space
ycheb = funeval(c,chebproj,x);      % values
dcheb = funeval(c,chebproj,x,1);    % first derivative
scheb = funeval(c,chebproj,x,2);    % second derivative

% Construct and evaluate cubic spline interpolant
spliproj = fundefn('spli',n,a,b);   % construct projection space
c     = funfitf(spliproj,'fapp06'); % project onto space
ycspl = funeval(c,spliproj,x);      % values
dcspl = funeval(c,spliproj,x,1);    % first derivative
scspl = funeval(c,spliproj,x,2);    % second derivative

% Plot function approximation error
figure(1);
subplot(2,1,1); plot(x,y-ycheb,'r'); ylabel('Chebychev');
subplot(2,1,2); plot(x,y-ycspl,'m'); ylabel('Cubic Spline');
subplot(2,1,1); title('Function Approximation Error');
subplot(2,1,2); xlabel('x');

% Plot first derivative approximation error
figure(2);
subplot(2,1,1); plot(x,d-dcheb,'r'); ylabel('Chebychev');
subplot(2,1,2); plot(x,d-dcspl,'m'); ylabel('Cubic Spline');
subplot(2,1,1); title('First Derivative Approximation Error');
subplot(2,1,2); xlabel('x');

% Plot second derivative approximation error
figure(3);
subplot(2,1,1); plot(x,s-scheb,'r'); ylabel('Chebychev');
subplot(2,1,2); plot(x,s-scspl,'m'); ylabel('Cubic Spline');
subplot(2,1,1); title('Second Derivative Approximation Error');
subplot(2,1,2); xlabel('x');

figure(4);
plot(x,y-ycheb)
title('Approximation Error');

prtfigs(mfilename,'Approximation Error for exp(-\alpha x)',4)
