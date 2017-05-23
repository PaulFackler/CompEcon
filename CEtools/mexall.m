% MEXALL Creates MEX files for CompEcon toolbox
% Run this function once after the toolbox is installed
% Requires a C compiler
% Type help mex for more information on creating mex files
% You may want to set compiler options to change speed/memory 
%   default settings
%
% This program need be run only once
%
% Note: this will compile all of the C files in the 
%   cetools and cetools\private directories

function mexall

% save the name of the current default directory
currentdir=cd;

% find and switch to the cetools directory
cetools=which('fundef');
while cetools(end)~='\' & cetools(end)~='/', cetools(end)=[];end
cd(cetools)

% mex all C files in the cetools directory
files=dir('*.c');
for i=1:length(files)
  %mex(files(i).name)
  eval(['mex -largeArrayDims ' files(i).name])
  disp(['mex file created for ' files(i).name])
end

% switch to the cetools\private directory
cd([cetools 'private'])
% mex all C files in the cetools\private directory
files=dir('*.c');
for i=1:length(files)
  %mex(files(i).name)
  eval(['mex -largeArrayDims ' files(i).name])
  disp(['mex file created for ' files(i).name])
end

% switch back to original default directory
cd(currentdir)