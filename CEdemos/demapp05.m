% DEMAPP05 Demonstrates alternative approximants for various functions
% Demonstrates Chebychev polynomial, cubic spline, and
% linear spline approximation for the following functions
%   1: y = sqrt(x+1);
%   2: y = exp(-x);
%   3: y = 1./(1+25*x.^2).
%   4: y = sqrt(abs(x));
close all

% Set endpoints of interpolation interval
a = -1;                       % left endpoint
b =  1;                       % right endpoint
n =  25;

% Construct Chebychev, cubic spline, and linear spline projection space
chebproj = fundefn('cheb',n,a,b);
csplproj = fundefn('spli',n,a,b);
lsplproj = fundefn('spli',n,a,b,1);

% Construct refined uniform grid for error ploting
x = nodeunif(2001,a,b);

for ifunc=1:4
   % Construct interpolants
   ccheb = funfitf(chebproj,'fapp05',ifunc);
   ccspl = funfitf(csplproj,'fapp05',ifunc);
   clspl = funfitf(lsplproj,'fapp05',ifunc);
   
   % Compute actual and fitted values on grid
   y = fapp05(x,ifunc);                    % actual
   ycheb = funeval(ccheb,chebproj,x);      % Chebychev
   ycspl = funeval(ccspl,csplproj,x);      % cubic spline
   ylspl = funeval(clspl,lsplproj,x);      % linear spline

   % Plot function approximation error
   figure(ifunc);
   ymin=floor(min(y)); ymax=ceil(max(y)); 
   if ifunc==3
     axislim=[a b -0.2 1.2];
   else
     axislim=[a b ymin ymax];
   end
   subplot(2,2,1); plot(x,y,'k');     axis(axislim);  title('Function');
   subplot(2,2,2); plot(x,ycheb,'k'); axis(axislim);  title('Chebychev');
   subplot(2,2,3); plot(x,ycspl,'k'); axis(axislim);  title('Cubic Spline');
   subplot(2,2,4); plot(x,ylspl,'k'); axis(axislim);  title('Linear Spline'); 
end

prtfigs(mfilename,'Approximation of $\exp(-x)$',2)
prtfigs(mfilename,'Approximation of $(1+25x^2)^{-1}$',3)
prtfigs(mfilename,'Approximation of $|x|^{0.5}$',4)