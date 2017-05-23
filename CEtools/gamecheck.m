% GAMECHECK Checks derivatives for dynamic game model function files
% USEAGE
%   err = gamecheck(model,s,x,e)
% INPUTS
%   model : a game model structure
%   s     : a matrix of state nodes
%   x     : a matrix of action values
%   e     : a matrix of shocks
% OUTPUT
%   err   : a 6x1 vector containing the maximal differences 
%           between explicit and numerical derivatives:
%             err1: df1/dx
%             err2: d^2f1/dx^2
%             err3: df2/dx
%             err4: d^2f2/dx^2
%             err5: dg/dx
%             err6: d^2g/dx^2
%
% See also: GAMESOLVE

% Copyright (c) 1997-2003, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function err=gamecheck(model,s,x,e)

[n,d]=size(s);
m=size(x,2);

if nargin<4; 
    if isfield(model,'w'), e=model.w'*model.e; 
    else,                  e=0;
    end
end

[f1,df1,d2f1]=feval(model.func,'f1',s,x,e,model.params{:});
[f2,df2,d2f2]=feval(model.func,'f2',s,x,e,model.params{:});
[g,dg,d2g]=feval(model.func,'g',s,x,e,model.params{:});

Jf1=fjac(model.func,3,'f1',s,x,e,model.params{:});
Jf2=fjac(model.func,3,'f2',s,x,e,model.params{:});
Jg=fjac(model.func,3,'g',s,x,e,model.params{:});

Hf1=fjac(model.func,[3 2],'f1',s,x,e,model.params{:});
Hf2=fjac(model.func,[3 2],'f2',s,x,e,model.params{:});
Hg=fjac(model.func,[3 2],'g',s,x,e,model.params{:});


ind1=find(~isinf(df1));
if isempty(ind1)
  err=NaN;
else
  err=max(abs(df1(ind1)-Jf1(ind1)));
end
ind2=find(~isinf(d2f1));
if isempty(ind2)
  err=[err NaN];
else
  err=[err max(abs(d2f1(ind2)-Hf1(ind2)))];
end

ind3=find(~isinf(df2));
if isempty(ind3)
  err=[err NaN];
else
  err=[err max(abs(df2(ind3)-Jf2(ind3)))];
end
ind4=find(~isinf(d2f2));
if isempty(ind4)
  err=[err NaN];
else
  err=[err max(abs(d2f2(ind4)-Hf2(ind4)))];
end

ind5=find(~isinf(dg));
if isempty(ind5)
  err=[err NaN];
else
  err=[err max(abs(dg(ind5)-Jg(ind5)))];
end

ind6=find(~isinf(d2g));
if isempty(ind6)
  err=[err NaN];
else
  err=[err max(abs(d2g(ind6)-Hg(ind6)))];
end

if length(ind1)+length(ind2)+length(ind3)+length(ind4)+length(ind5)+length(ind6) ...
    <6*size(s,1)
  warning('Infinite derivatives encountered')
end

if max(err)>1.e-4
   disp('Possible Error in Derivatives')
   disp('Discrepancies in derivatives = ')
   disp(err)
end