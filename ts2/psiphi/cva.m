function [A,scales,yp,yt,ep,pl] = ...
    cva(z,np,ps);

% [A,sc,yp,yt,ep,pl]=cva(data,Np,ps)  cva interpolation for variable step size model.
%
% The variable step size model F(x(t);tau)=x(t+tau)+noise (from continous)
% is interpolated. That is; let y be the next k observations and x the
% current one. z=(F(x(t),tau(1)),...F(x(t),tau(??))) are all the predictions
% of future observations from a variable step size model. Then y may be
% estimated as follows;
%     y=Az+e
%     z=(F(x(t),tau(1)),...F(x(t),tau(??)))
% x are the time delay embedding of data.
% Instead of taking just (F(x(t),tau(1)),...F(x(t),tau(??))) one may take
% polynomial terms with these factors.
% Np is the (maxmimum) order of the polynomials.
% ps is the maximum prediction step (i.e. k).
% sc are the normalising factors associated with the calculation of A.
%
% I'm not quite sure what is going on, and I certainly don't understand the
% theory. It is somehow related to Canonical Variate Analysis (and therefore
% SVD) and appears to provide a linear, quadratic or spline interpolation of 
% (F(x(t),tau(1)),...F(x(t),tau(??))) in the case when diff(tau) are large.
%

% Copyright (c) 1998 by Michael Small.
%
% Please see the copyright notice included in this distribution
% for full details.
%
%
% File   cva.m
%   $Id$
%
% Created by Michael Small (<watchman@>) on Fri Aug  7 1998
%
% $Log$

if nargin<3
  ps=[];
end;
if nargin<2,
  np=1;
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

disp('fitting');
%calculate A matrix
[x,scales]=normalize(x);
resc=RMS(y(1,:))/RMS(x(1,:));
scales=scales*resc;
x=x/resc;
x=x';
A=y/x;

disp('testing');
%do some testing
yp=A*x;
yt=y;
ep=yp-yt;
% End of cva.m
