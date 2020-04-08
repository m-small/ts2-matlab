function [A,scales,yp,yt,ep,pl] = ...
    cva2(x,y,np,nf,ps);

% [A,sc,yp,yt,ep,pl]=cvaf(z,y,Np,ps)  cva interpolation for variable step size model.
%
% same as cva, but y and z=F(x) are specified.
% y and z are both matrices each row of z (column of y) corresponds to one
% embedded point (future vector). They are back-the-front from what you get
% out of embed. i.e. higher ordinate matrix elements are further into the
% future. Its ass-over.
%
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

if nargin<4,
  ps=[];
end;
if nargin<3,
  nf=100;
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

disp('multiplying');
%predict
pl=x;

%powers of x
[x,pie]=polyn(x,np);

%generate the future vectors/future matrix
F=fseries(ps+1,nf);
F(ps+1,:)=[];

disp('fitting');
%calculate A matrix
[x,scales]=normalize(x);
resc=RMS(y(1,:))/RMS(x(1,:));
scales=scales*resc;
x=x/resc;
x=x';
mf=mean(diag(F*F'));
A=(F'*y/mf)/x;

disp('testing');
%do some testing
yp=(F*A)*x;
yt=y;
ep=yp-yt;
% End of cva.m










