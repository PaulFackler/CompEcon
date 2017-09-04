%% PRINTFIGURES
%
% Save figures for overheads and manuscripts

%  Copyright(c) 1997-2016
%   Mario J. Miranda - miranda.4@osu.edu
%   Paul L. Fackler  - paul_fackler@ncsu.edu

function printfigures(nam,titon)

return

h = findobj('type','figure');
for i=1:length(h)
  figure(i)
  if nargin<2||~titon==1
    title([])
  end
  box off
  eval(['print -depsc ' ['figures\' nam num2str(i)]])
end