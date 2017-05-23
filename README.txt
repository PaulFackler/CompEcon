Revised: 1/9/2011

If you are reading this file, you must have unzipped the file
compecon.zip. This should be unzipped into a directory called 
COMPECON, which itself has two subdirectories
CETOOLS and CEDEMOS (directory information attached to the files) 

To use the toolbox, you must change the MATLAB path so it can
find the files. Type 
  cepath='d:\compecon\'; path([cepath 'cetools;' cepath 'cedemos'],path);
at the MATLAB command line. It is assumed that COMPECON is on the D: drive;
if not, change cepath in the line above. To have this occur automatically add 
this line to your matlabrc.m or startup.m files (these can be located using the which command). 

It is preferable to locate the compecon directory as a subdirectory of main 
matlab path; this leads MATLAB to skip checking file time stamps and leads to speedier execution.

To see a listing of all the toolbox functions, type help cetools. 
To see a listing of all the demonstration files, type help cedemos. 

This toolbox is available on the internet at
  www4.ncsu.edu/~pfackler/compecon
Please do not distribute this toolbox; instead provide the URL.

IMPORTANT NOTE on MEX files
A number of the routines provided in this toolbox are coded in C
as MEX files. We have provided versions of these for Windows users of
MATLAB ver. 7. These may not work on your system however. Most but not all
of the functionality of the toolbox will run without MEX files but
a number of procedures will run slowly and may encounter memory limitations
using the m-file versions.

The utility MEXALL will create these mex files.
It is assumed that you have a C compiler and can create mex files. See
MATLAB's Applications Interface documentation for details. 32-versions of Matlab 
come supplied with the LCC compilers. The first time you run MEXALL you may be asked
to chose a compiler. 64-bit versions do not ship with a compiler and the MEX files are not
tested for 64-bit platforms in any case. Someone willing to test these on a 64-bit platform 
should contact me at paul_fackler@ncsu.edu.