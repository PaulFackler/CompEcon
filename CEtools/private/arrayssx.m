% ARRAYSSX Utility used with ARRAYSS
% USAGE
%    [fxnew,ff,aa]=arrayssx(x,a,b,fx);  
% INPUTS
%    x : an evaluation point
%    a : lower bounds
%    b : upper bounds
%   fx : function values at x
% OUTPUTS
%   fxnew : value of semi-smooth function at x
%      ff : derivative of function at x
%      aa : derivative of function at x
%
% The reformulation uses
%   phi^-(phi^+(fx,a-x),b-x)
% where
%   phi^+(y,z)=y+z+sqrt(y.^2+z^2)
%   phi^-(y,z)=y+z-sqrt(y.^2+z^2)


function [fxnew,ffout,aaout]=arrayssx(x,a,b,fx)
  n=numel(x);
  fxnew=zeros(size(x));
  if nargout>1, ffout=zeros(size(x));
    if nargout>2, aaout=zeros(size(x)); end
  end  
  if length(a)==1, a=a+zeros(size(x)); end
  if length(b)==1, b=b+zeros(size(x)); end
  for j = 1:n
    % compute phi+ 
    if isinf(a(j)), d=fx(j); 
    else
      da=a(j)-x(j);
       if abs(f(j))>abs(da), y=fx(j); z=da;
       else                  y=da;    z=fx(j);
       end
       if y==0, d=0;
       else
         z=z/y;
         dplus=sqrt(1+z*z);
         if y>0, d=y*(1+dplus+z);
         else    d=y*(z-((1-dplus)*(1-dplus)+z*z)/dplus/2);
         end
       end
    end
    % compute phi- 
    if isinf(b(j)), fxnew(j)=d;
    else
      db=b(j)-x(j);
      if abs(d)>abs(db), g=d;  h=db;
      else               g=db; h=d; 
      end
      if g==0, fxnew(j)=0; 
      else
        h=h/g;
        dminus=sqrt(1+h*h);
        if (g<0), fxnew(j)=g*(1+dminus+h);
        else      fxnew(j)=g*(h-((1-dminus)*(1-dminus)+h*h)/dminus/2);
        end
      end
    end
    % compute Jacobian factors if requested 
    if nargout>1
      if isinf(b(j))
        ff=1; aa=1; bb=0;
      else  
        if g<0, dminus=-dminus; end
        temp1=1-1/dminus;  temp2=1-h/dminus;
        if abs(d)>abs(db), ff=temp1; aa=temp1; bb=temp2;
        else              ff=temp2; aa=temp2; bb=temp1;
        end
      end
      if isinf(a(j)), aa=0;
      else
        if y<0, dplus=-dplus; end
        temp1=1+1/dplus; temp2=1+z/dplus;
        if abs(f(j))>abs(da), ff = ff*temp1; aa=aa*temp2;
        else                  ff = ff*temp2; aa=aa*temp1;
        end
      end
      ffout(j)=ff;
      if nargout>2, aaout(j)=aa+bb; end
    end
  end