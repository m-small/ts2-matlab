function base=rb_stretch_all(base,X,v,func,e,meth);

% function [base,phi,scale,imprv]=rb_stretch(nbase,X,v,func,e,meth);
%
% perform Nelder-Mead Simplex minimization to minimize -abs(phi'*e) subject to the params
% of each single basis function in nbase. This uses MATLABs fminsearch (was fmins in v4).
% X,v,func,e, and meth all come from rb_topdown
%
% NB: we have
%  phi=rb_Phi(X,base,v,func);
%  [phi,scale]=normalize(phi);
%  obj=-abs(phi'*e);
% obj is the function to be minimized.
%
% M. Small 
% Created: 6/3/02
% Updated: 6/3/02

alpha=0.1;

opt=optimset(optimset,'Display','off','TolX',1e-3); 

[d,nc]=size(base.centres);

for i=1:nc,
  vv= v(base.strategy(i),:);
  is= find(vv>=0);
  vs= 1+vv(is);

  x0=[base.centres(vs,i);base.radii(:,i)];

  ind=find(meth=='c' | meth=='l');
  meth(ind)=[];

  xi=fminsearch('rb_wobbly',x0,opt,X,vv,func(base.func),e,meth,alpha);
  imprv=all(xi==x0); %if imrpv then no imprv!

  lvs=length(vs);
  base.centres(vs,i)=xi(1:lvs);
  base.radii(:,i)=xi((lvs+1):end);

end;

phi=rb_Phi(X,base,v,func,meth);
phi=phi(:,end);
[phi,scale]=normalize(phi);
scale=scale;
