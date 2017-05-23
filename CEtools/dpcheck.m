% DPCHECK Checks derivatives for discrete time DP model function files
% USEAGE
%   err = dpcheck(model,s,x,e)
% INPUTS
%   model : a DP model structure
%   s     : a matrix of state nodes
%   x     : a matrix of action values
%   e     : a matrix of shocks
% OUTPUT
%   err   : a 4x1 vector containing the maximal differences 
%           between explicit and numerical derivatives:
%             err1: df/dx
%             err2: d^2f/dx^2
%             err3: dg/dx
%             err4: d^2g/dx^2
%
% See also: DPSOLVE

% Copyright (c) 1997-2002, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function err=dpcheck(model,s,x,e)

[n,d]=size(s);
m=size(x,2);

if nargin<4; 
    if isfield(model,'w'), e=model.w'*model.e; 
    else                   e=0;
    end
end

Jf=fjac(model.func,3,'f',s,x,e,model.params{:});
%if size(Jf,1)>1 || size(Jf,2)>1, Jf=diag(Jf); end
Jg=fjac(model.func,3,'g',s,x,e,model.params{:});
%if size(Jg,1)>1 || size(Jg,2)>1, Jg=diag(Jg); end

Hf=fjac(model.func,[3 2],'f',s,x,e,model.params{:});
Hg=fjac(model.func,[3 2],'g',s,x,e,model.params{:});

[f,df,d2f]=feval(model.func,'f',s,x,e,model.params{:});
[g,dg,d2g]=feval(model.func,'g',s,x,e,model.params{:});

ind1=find(~isinf(df));
if isempty(ind1)
  err=NaN;
else
  err=max(abs(vec(df(ind1))-vec(Jf(ind1))));
end
ind2=find(~isinf(d2f));
if isempty(ind2)
  err=[err NaN];
else
  err=[err max(abs(vec(d2f(ind2))-vec(Hf(ind2))))];
end
ind3=find(~isinf(dg));
if isempty(ind3)
  err=[err NaN];
else
  err=[err max(abs(vec(dg(ind3))-vec(Jg(ind3))))];
end
ind4=find(~isinf(d2g));
if isempty(ind4)
  err=[err NaN];
else
  err=[err max(abs(vec(d2g(ind4))-vec(Hg(ind4))))];
end


if 0 % commented out because derivatives with respect to s are no longer used
% Check derivatives with respect to s
h = eps^(1/3)*max(abs(s),1);
sh1=s+h; sh0=s-h;
h=sh1-sh0;
Jf=zeros(n,d);
Jg=zeros(n,d,d);
for j=1:d;
   ss = s; 
   ss(:,j) = sh1(:,j);  
   f1=feval(model.func,'f',ss,x,e,model.params{:});
   g1=feval(model.func,'g',ss,x,e,model.params{:});
   ss(:,j) = sh0(:,j); 
   f0=feval(model.func,'f',ss,x,e,model.params{:});
   g0=feval(model.func,'g',ss,x,e,model.params{:});
   Jf(:,j) = (f1-f0)./h(:,j);
   Jg(:,:,j) = (g1-g0)./h(:,j+zeros(1,d));
end

ind5=find(~isinf(dsf));
if isempty(ind5)
  err=[err NaN];
else
  err=[err max(abs(vec(dsf(ind5))-vec(Jf(ind5))))];
end
ind6=find(~isinf(dsg));
if isempty(ind6)
  err=[err NaN];
else
  err=[err max(abs(vec(dsg(ind6))-vec(Jg(ind6))))];
end
end

if length(ind1)+length(ind2)+length(ind3)+length(ind4) ...
    <4*size(s,1)
  warning('Infinite derivatives encountered')
end

if max(err)>1.e-4
   disp('Possible Error in Derivatives')
   disp('Discrepancies in derivatives = ')
   disp(err)
end

function v=vec(A)
v=A(:);