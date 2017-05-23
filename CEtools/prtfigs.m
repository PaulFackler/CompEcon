% PRTFIGS Utility to save figures as .eps files
% USAGE
%   prtfigs(nam,TitleSize,LabelSize)
% To activate use:
%   optset('prtfigs','dir','d:\booknew\figs\')
% (change the path name as appropriate)

% Copyright (c) 1997-2000, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function prtfigs(nam,TitleText,figlist)

global prtfigs_options

if isempty(prtfigs_options) | isempty(prtfigs_options.dir)
  return
end

%savefigs(nam)

warnings=optget('prtfigs','warnings',0);

%savehandle = findobj('buttondownfcn','moveaxis');
%set(savehandle,'buttondownfcn','')




% PRTFIGS Prints the currently open figures to files
% names NAM1, NAM2, etc., where NAM is passed to PRTFIGS


%set(0,'DefaultAxesPosition',[0.1    0.25    0.55    0.6667]);
%set(0,'DefaultFigurePaperPosition',[1.5 3.625 5.5 3.85]);
%set(0,'DefaultAxesFontName','Times-Roman');
%set(0,'DefaultAxesFontSize',12);


if nargin<3 
  disp('Not enough inputs in prtfigs')
  return 
end
if ~exist('figlist','var')
  disp('No figures printed')
  return
end

LabelSize=20;
TitleSize=20;
LabelFont='Times-Roman';

hassubplots=0;
for i=1:length(figlist);

  h=findobj(figlist(i),'type','axes');
  
  nplots=length(h);

  for i=length(h):-1:1; 
    if strcmp(get(h(i),'Tag'),'legend')
      nplots=nplots-1;
      set(findobj(h(i),'type','text'),'FontName',LabelFont);
    else
      set(get(h(i),'XLabel'),'FontSize',LabelSize,'FontName',LabelFont);
      set(get(h(i),'YLabel'),'FontSize',LabelSize,'FontName',LabelFont);
      set(get(h(i),'ZLabel'),'FontSize',LabelSize,'FontName',LabelFont);
    end
  end

  if length(figlist)>1 & nplots>1
    error('plots with subplots must be saved separately')
  end

  if nplots>1, hassubplots=1; end

end

if hassubplots & length(h)<=4
    h=sort(h);
    n=num2str(figlist(1));
    if length(h)==2
      texline=['\figsubtwo{' nam '}{' n '}'];
    else  
      texline=['\figsubs{' nam '}{' n '}'];
    end
    texnames=[];
    nam=[prtfigs_options.dir nam];
    texline=[texline '{' TitleText '}'];
    figure(figlist(1));
    for j=1:length(h)
      if ~strcmp(get(h(i),'Tag'),'legend')
        th=get(h(j),'title');
        texline=[texline '{' get(th,'string') '}'];
        axes(h(j));
        th=title(char(96+j));
        set(th,'fontname',LabelFont,'fontsize',LabelSize/2)
      end
      set(get(h(j),'XLabel'),'FontSize',LabelSize/2,'FontName',LabelFont);
      set(get(h(j),'YLabel'),'FontSize',LabelSize/2,'FontName',LabelFont);
      set(get(h(j),'ZLabel'),'FontSize',LabelSize/2,'FontName',LabelFont);
      %tp=get(th,'position');
      %te=get(th,'extent');
      %tp(2)=tp(2)-te(4)/2;
      %set(th,'position',tp);
      set(h(j),'Fontsize',get(h(j),'FontSize')/2);
      set(findobj(gca,'type','text'),'FontSize',get(0,'DefaultTextFontSize')/2);
    end
    eval(['print -f' n ' -depsc  -cmyk -loose ' nam n])
    disp(texline)
  
else
  switch length(figlist)
    case 1
    texline=['\figone{' nam '}'];
    case 2
    texline=['\figtwo{' nam '}'];
    case 3
    texline=['\figthree{' nam '}'];
    case 4
    texline=['\figfour{' nam '}'];
  end

  texnames=[];
  nam=[prtfigs_options.dir nam];
  
  if length(figlist)>1
    for i=1:length(figlist)
      figure(figlist(i))
      texline=[texline '{' num2str(figlist(i)) '}'];
      if i==1,
        texnames=[char(96+i) '. ' get(get(gca,'Title'),'String')];
      elseif i==length(figlist) 
        if rem(i,2)
          texnames=[texnames ' \\ \multicolumn{2}{c}{' char(96+i) '. ' get(get(gca,'Title'),'String') '}'];
        else
          texnames=[texnames ' & ' char(96+i) '. ' get(get(gca,'Title'),'String')];
        end
      else
        if rem(i,2)
          texnames=[texnames ' \\ ']; 
        else
          texnames=[texnames ' & ']; 
        end
        texnames=[texnames  char(96+i) '. ' get(get(gca,'Title'),'String')];
      end
      h=title(char(96+i));
      set(h,'FontName',LabelFont,'FontSize',TitleSize) 
      n=num2str(figlist(i));
      eval(['print -f' n ' -depsc   -cmyk ' nam n])
    end
    texline=[texline '{' TitleText '}{' texnames '}'];
  else
    figure(figlist(1))
    n=num2str(figlist(1));
    texline=[texline '{' n '}{' TitleText '}'];
    title('');
    eval(['print -f' n ' -depsc  -cmyk -loose ' nam n])
  end
  disp(texline)
end 



