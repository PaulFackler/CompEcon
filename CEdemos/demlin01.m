% DEMLIN01 Compares effort to solve Ax=b by different methods
close all

disp(' ')
disp('DEMLIN01 Compares effort to solve Ax=b by different methods')

disp(' ')
disp('Approximate flops required to solve n by n linear')
disp('equation Ax=b m times using A\b and')
disp('inv(A)*b, computing inverse only once.')
disp(' ')
disp('    m       n        A\b     A^-1*b')
disp(' ')

for m=[1,100]
for n=[50,500]

   A = rand(n,n);
   b = rand(n,1);

   tic
   for j=1:m
      x = A\b;
   end
   f1 = toc;

   tic
   Ainv = inv(A);
   for j=1:m
      x = Ainv*b;
   end
   f2 = toc;

   fprintf('  %3i     %3i  %9.2f  %9.2f \n',m,n,f1,f2)

end
end