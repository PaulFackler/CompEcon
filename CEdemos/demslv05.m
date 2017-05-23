% DEMSLV05 Compare various linear complementarity algorithms
disp('DEMSLV05 Compare various linear complementarity algorithms')
close all

% Generate test problem data
n  = 8;
z = randn(n,2)-1;
a = min(z,[],2);
b = max(z,[],2);
q  = randn(n,1);
M  = randn(n,n);
M  = -M'*M;
xinit  = randn(n,1);

warning off;

fprintf('Solution to Randomly Generated LCP\n');
fprintf('name                    time      norm        x1       x2\n');
for i=1:4
   tic
   switch i
   case 1
      name = 'Newton semismooth     ';
      optget('lcpsolve','type','smooth');
      [x,z] = lcpsolve(M,q,a,b,xinit);
   case 2
      name = 'Newton minimax        ';
      optget('lcpsolve','type','minmax');
      [x,z] = lcpsolve(M,q,a,b,xinit);
   case 3
      name = 'Baard                 ';
      [x,z] = lcpbaard(M,q,a,b,xinit);
   case 4
      name = 'Lemke                 ';
      [x,z] = lcplemke(M,q,a,b,xinit);
   end 
   err = norm(minmax(x,a,b,z));
   fprintf('%s%6.2f  %8.1e  %8.3f %8.3f\n',name,[toc err x(1) x(2)])
end

warning on;
