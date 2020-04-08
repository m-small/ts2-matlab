% =setpath()  set matlab path to use cva programs and models
%
% it works in UNIX.

% Copyright (c) 1998 by Michael Small.
%
% Please see the copyright notice included in this distribution
% for full details.
%
%
% File   setpath.m
%   $Id$
%
% Created by Michael Small (<watchman@>) on Thu Sep 10 1998
%
% $Log$

home=getenv('HOME');
p=[home,'/matlab/cva/:',home,'/matlab/cva/models/:',home,'/matlab/cva/other/'];
path(p,path);


% End of setpath.m
