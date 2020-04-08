function base=rb_stretch_r(base,X,v,func,e,meth);

% function base=rb_stretch(nbase,X,v,func,e,meth);
%
% perform Nelder-Mead Simplex minimization to minimize -abs(phi'*e) subject to the radii
% of each singla basis function in nbase. This uses MATLABs fminsearch (was fmins in v4).
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



opt=optimset(optimset,'Display','off','TolX',1e-3); 

[d,nc]=size(base.centres);

for i=1:nc,
  vv= v(base.strategy(i),:);
  is= find(vv>=0);
  vs= 1+vv(is);

  xp=base.centres(vs,i);
  x0=base.radii(:,i);

  ind=find(meth=='c' | meth=='l');
  meth(ind)=[];

  xi=fminsearch('rb_stretchy',x0,opt,xp,X,vv,func(base.func),e,meth);

  base.radii(:,i)=xi;

end;

