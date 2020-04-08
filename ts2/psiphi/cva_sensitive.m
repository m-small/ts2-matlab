function [A,scales,prc,yp,yt,ep,pl] = ...
    cva_sensitive(x,y,np,ps);

% [A,sc,prc,yp,yt,ep,pl]=cva2(z,y,Np,ps)  cva interpolation for variable step size model.
%
% same as cva, but y and z=F(x) are specified.
% y and z are both matrices each row of z (column of y) corresponds to one
% embedded point (future vector). They are back-the-front from what you get
% out of embed. i.e. higher ordinate matrix elements are further into the
% future. Its ass-over.
%
% Furthermore sensitive is employed to calculated the MDL best A (row-wise!)
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
% prc is the matrix of precisions corresponding (elemnetwise) to A.
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

disp('multiplying');
%predict
pl=x;

%powers of x
[x,pie]=polyn(x,np);

disp('fitting');
%calculate A matrix
[x,scales]=normalize(x);
resc=RMS(y(1,:))/RMS(x(1,:));
scales=scales*resc;
x=x/resc;
x=x';
%do A=y/x; with sensitive
[nx,lx]=size(x);
A=zeros(nx,ps);prc=zeros(nx,ps);
[level,is_on]=trace;
trace('off');
timer('s');
for i=1:ps
  timer(i/ps);
  [best_basis,lambda,precision]= sensitive(x',y(i,:));
  prc(best_basis,i)=precision;
  A(best_basis,i)=lambda;
end;timer('f');
trace(is_on,level);
A=A';prc=prc';

disp('testing');
%do some testing
yp=A*x;
yt=y;
ep=yp-yt;
% End of cva.m










