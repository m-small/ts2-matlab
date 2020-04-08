function [yp,yt,ep,pl] = ...
    cva_pred(z,A,sc,np,ps);

% [yp,yt,ep,pl] = ...
%    cva_pred(z,A,sc,np,ps);
%
% predictions on z, using cva matrix A (from cva.m) to
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

if na<5
  ps=[];
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


disp('predicting');
%predict
y=[];x=[];dw=[0 cumsum(-diff(w))];dw=fliplr(dw);
for i=1:length(w),
  [yp,yt]=dc_predict(z,w(i));
  x=[x yp(1:(-dw(i)+length(yp)))];
end;
lx=length(x(1,:));
x=x(1:(length(yp)-(ps-maxw)),:);
pl=x;

disp('building');
%how long should the prediction thing be?
maxv=max(v(~isnan(v(:))))-1;
lp=length(z)-maxv-ps;
for wi=min(-w):ps,%max(-w),
  y(wi+1,:)=z(maxv+wi+(1:lp));
end;

%powers of x
[x,pie]=polyn(x,np);
x=x./(ones(length(x),1)*sc);
disp('testing');
%do some testing
yp=A*x';
yt=y;
ep=yp-yt;

% End of cva_pred.m
