function [yp,pl] = ...
    cva_predf(z,A,sc,np,nf,ps);

% [yp,pl] = ...
%    cva_predf(y,A,sc,np,nf,ps);
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

if na<6
  ps=[];
end;
if na<5
  nf=100;
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

if isempty(ps),
  ps=maxw;
end;


%powers of x
[x,pie]=polyn(z,np);
x=x./(ones(length(x),1)*sc);

%generate the future vectors/future matrix
F=fseries(ps+1,nf);
F(ps+1,:)=[];

disp('testing');

%do some testing
yp=(F*A)*x';
% End of cva_pred.m
