% Makefile.
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 03-20-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 05-05-2013

path0 = cd;

cd 'lib/cell';
mex GCC='/usr/bin/gcc-6' cellss.cpp;
mex GCC='/usr/bin/gcc-6' oness.cpp;
mex GCC='/usr/bin/gcc-6' zeross.cpp;
cd(path0);

cd 'lib/text';
mex GCC='/usr/bin/gcc-6' atoi.cpp;
mex GCC='/usr/bin/gcc-6' atof.cpp;
mex GCC='/usr/bin/gcc-6' tokenise.cpp;
cd(path0);

cd 'lib/img';
mex GCC='/usr/bin/gcc-6' maskOver.cpp;
cd(path0);

cd 'src/ali/dtw';
mex GCC='/usr/bin/gcc-6' dtwFord.cpp;
mex GCC='/usr/bin/gcc-6' dtwBack.cpp;
mex GCC='/usr/bin/gcc-6' dtwFordAsy.cpp;
cd(path0);

cd 'src/ali/help';
mex GCC='/usr/bin/gcc-6' rowBd.cpp;
cd(path0);

cd 'src/ali/imw';
mex GCC='/usr/bin/gcc-6' timewarp.cpp;
cd(path0);
