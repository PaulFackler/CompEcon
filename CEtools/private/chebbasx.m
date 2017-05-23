% CHEBBASX A utility used by CHEBBAS
% Coded as a C-mex file for additional speed

% Copyright (c) 1997-2001, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu
function bas=chebbasx(n,a,b,x)
 
global CompEcon_MEXwarned
 
if isempty(CompEcon_MEXwarned)
  disp('Warning: You are using the m-file version of a function that is coded as a C Mex file.')
  disp('  Running this M file version may be significantly slower and more memory intensive.')
  disp('  Use MEXALL to create the executable (MEX or DLL) and make sure it is on the MATLAB path.')
  CompEcon_MEXwarned=1;
end

  z = (2/(b-a))*(x-(a+b)/2);
  m=size(z,1);
  bas=zeros(m,n);
  bas(:,1)=ones(m,1);
  bas(:,2)=z;
  z=2*z;
  for i=3:n
    bas(:,i)=z.*bas(:,i-1)-bas(:,i-2);
  end
