function [mu]=rb_stretchy(xi,xc,X,v,func,e,meth);

% function [mu,phi]=rb_stretchy(xi,xc,X,v,func,e,meth);
%
% evaluation of obj=mu=-phi'*e
%
% see rb_wobble
%
% modified to incorporate rh_wobbly_Phi
%
% M. Small 
% Created: 26/6/99
% Updated: 11/2/02

beta=9;

is= find(v>=0);
vs= 1+v(is);

base.centres(vs,1)=xc;
base.radii=xi;
base.strategy=1;
base.func=1;
phi=rb_Phi(X,base,v,func,meth);

phi=normalize(phi);
mu=-abs(phi'*e); %this is the sensitivity





