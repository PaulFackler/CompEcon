% DEMLIN02 Demonstrates ill-conditionsing of Vandermonde matrix
close all

disp(' ')
disp('DEMLIN02 Demonstrates ill-conditionsing of Vandermonde matrix')

disp('     Conditioning of nxn')
disp('     Vandermonde matrix')
disp(' ')
disp('     Norm of    Condition  ')
disp(' n   I-V\V      Number of V')
disp(' ')

warning off
for n=5:5:50
   v = vander([1:n]'); errv = norm(eye(n,n)-v\v); conv = cond(v);
   fprintf ('%2i %10.1e %10.1e\n',n,errv,conv);
end
warning on
