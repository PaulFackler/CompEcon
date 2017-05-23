% DEMSLV06 Demonstrates different NCP methods
function demslv06
close all

disp(' ')
disp('DEMSLV06 Demonstrates different NCP methods')

% Generate problem test data
z = rand(2,2)*10;
a = min(z,[],2);
b = max(z,[],2);
xinit = (b+a)/2;
xinit=rand(2,1);

warning off;

fprintf('Solution to Randomly Generated NCP\n');
fprintf('name                      time   norm     x(1)  x(2)\n');
for i=1:6
   tic
   switch i
   case 1
      name = 'Newton semismooth     ';
      optset('ncpsolve','type','smooth');
      %optset('ncpsolve','maxsteps',0);
      [x,z] = ncpsolve('fslv06',a,b,xinit);
   case 2
      name = 'Newton minimax        ';
      optset('ncpsolve','type','minmax');
      %optset('ncpsolve','maxsteps',0);
      [x,z] = ncpsolve('fslv06',a,b,xinit);
   case 3
      name = 'SLCP Newton semismooth';
      [x,z] = ncpjose('smooth','fslv06',a,b,xinit);
   case 4
      name = 'SLCP Newton minimax   ';
      [x,z] = ncpjose('minmax','fslv06',a,b,xinit);
   case 5
      name = 'SLCP Baard            ';
      [x,z] = ncpjose('baard','fslv06',a,b,xinit);
   case 6
      name = 'SLCP Lemke            ';
      [x,z] = ncpjose('lemke','fslv06',a,b,xinit);
   end 
   err = norm(minmax(x,a,b,z),inf);
   fprintf('%s%8.2f%10.1e%6.2f%6.2f\n',name,[toc err x(1) x(2)])
end

warning on;

optset('ncpsolve','defaults');