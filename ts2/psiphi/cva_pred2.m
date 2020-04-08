function [yp,pl] = ...
    cva_pred(z,A,sc,np,xx);

% [yp,pl] = ...
%    cva_pred2(y,A,sc,np,x);
%
% predictions on y==y_pred, using cva matrix A (from cva.m) to
% smooth/interpolate/extrapolate with np-pth order polynomials up ps points
% in the future. Typically one would do something like
%   continuous(z,v,w,...);
%   [A,yp,yt,ep]=cva(z,np,ps);
%   [yp,yt,ep,pl]=cva_pred(z_test,A,np,ps];
% (where, probably, w is 1/4 period predictions ahead)
%
% yp -predictions, yt - test values, ep - error (yp-yt).
%
% if x is given the A is fit of yp=A[polyn(y) x];

% Copyright (c) 1998 by Michael Small.
%
% Please see the copyright notice included in this distribution
% for full details.
%
%
% File   cva_pred.m
%   $Id$
%
% Created by Michael Small (<watchman@>) on Mon Aug 10 1998
%
% $Log$

na=nargin;

if na<5,
  xx=[];  
end;
if na<4,
  np=1;
end;
if na<3,
  sc=1;
end;

%get globals
pl_globals

v=pl_embed;
w=pl_pred_vect;
dist=pl_timescale;
vmax=max(v(~isnan(v)));
maxw=-min(w);



%powers of x
[x,pie]=polyn(z,np);
x=[x xx'];
x=x./(ones(length(x(:,1)),1)*sc);
disp('testing');
%do some testing
yp=A*x';
% End of cva_pred.m
