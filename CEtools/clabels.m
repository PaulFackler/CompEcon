% CLABELS Returns labels from a contour plot matrix
% USAGE
%   [labels,x,y]=clabels(c);
% Also returned are the values of (x,y) at the last point plotted
% c is a contour matrix as described in documentation on the CONTOURC
% function. It is obtained using
%   c=contour(x,y,z);

% Copyright (c) 1997-2000, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function [labels,x,y]=clabels(c)

n=size(c,2);
ind=1;
i=1;
while ind<n
  z=c(1,ind);
  ind=ind+c(2,ind)+1;
  x(i)=c(1,ind-1);
  y(i)=c(2,ind-1);
  temp=num2str(z);
  labels(i,1:length(temp))=temp;
  i=i+1;
end
x=x(:);
y=y(:);