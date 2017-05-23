% DEMAPP03 Compare conditioning of Vandermonde and Chebychev matrices
disp('DEMAPP03 Compare conditioning of Vandermonde and Chebychev matrices')
close all
warning off

disp(' ')
disp('        Conditioning of Vandermonde matrix (V) and       ')
disp('             Chebychev matrix (B) of order n             ')
disp(' ')
disp('     Norm of    Condition      Norm of    Condition      ')
disp(' n   I-V\V      Number of V    I-B\B      Number of B    ')
disp(' ')

for n=5:5:50
   V = vander([1:n]');
   S = fundefn('cheb',n,0,1);
   B = funbas(S);
   errv = norm(eye(n,n)-V\V);
   errc = norm(eye(n,n)-B\B);
   fprintf ('%2i %10.1e %10.1e     %10.1e %10.1e \n', ...
            n,errv,cond(V),errc,cond(B));
end
warning on
