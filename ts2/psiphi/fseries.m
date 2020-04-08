function F = ...
    fseries(n,ord)

% F=fseries(n,ord)  n point ord order Fourier series
%
% n point ord order sine/cosine series Fourier series
% F has the form
% F=[ones sin(t) cos(t) sin(2t) cos(2t) .... sin(ord.t) cos(ord.t)];


% Copyright (c) 1998 by Michael Small.
%
% Please see the copyright notice included in this distribution
% for full details.
%
%
% File   fseries.m
%   $Id$
%
% Created by Michael Small (<watchman@>) on Wed Sep  2 1998
%
% $Log$

F=ones(1,n);
ft=0:(2*pi/(n-1)):(2*pi);

for i=1:ord,
  F=[F; sin(ft*i); cos(ft*i)];
end;

F=F';


% End of fseries.m
