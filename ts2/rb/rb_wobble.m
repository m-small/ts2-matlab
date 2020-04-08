function [base,phi,scale,imprv]=rb_wobble(base,X,v,func,e,meth,alpha);

% function [base,phi,scale,imprv]=rb_wobble(base,X,v,func,e,meth,alpha);
%
% perform Nelder-Mead Simplex minimization to minimize -abs(phi'*e) subject to the parameters
% of base, a single basis function. This uses MATLABs fminsearch (was fmins in v4).
% X,v,func,e, and meth all come from rb_topdown
% alpha is the penalty to avoid getting centres too far from the
% data (alpha>>0 means centres may be far from the data, try alpha=0.1)
%
% NB: we have
%  phi=rb_Phi(X,base,v,func);
%  [phi,scale]=normalize(phi);
%  obj=-abs(phi'*e);
% obj is the function to be minimized.
%
% imprv is 1 if wobbling helped, 0 otherwise.
%
% M. Small 
% Created: 26/6/99
% Updated: 11/2/02

if nargin<7,
  alpha=0.1;
end;

opt=optimset(optimset,'Display','off','TolX',1e-3); 
% set a low tolerance to prevent "overfitting"

vv= v(base.strategy,:);
is= find(vv>=0);
vs= 1+vv(is);

x0=[base.centres(vs);base.radii;];

ind=find(meth=='c' | meth=='l');
meth(ind)=[];

xi=fminsearch('rb_wobbly',x0,opt,X,vv,func(base.func),e,meth,alpha);
imprv=all(xi==x0); %if imrpv then no imprv!

lvs=length(vs);
base.centres(vs)=xi(1:lvs);
base.radii=xi((lvs+1):end);

phi=rb_Phi(X,base,v,func,meth);
phi=phi(:,end);
[phi,scale]=normalize(phi);
scale=scale;


