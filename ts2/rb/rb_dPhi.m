function T = rb_dPhi(x,base,v,func)

% X = rb_dPhi(x,base,v,func)
% 
% make an array X  at centres base.centre of radius base.radii 
% using embedding strategy base.strategy of v and basis function base.func.
% X(i,:,:) is the deriviative of of X=rb_Phi(x,base,v,func) with respect 
% to the i-th co-ordinate.
%
% X is dx-by-nx-by-(dx+1+nc) where  [dx,nx]=size(x) and nc= # of r.b. func.
%
% M. Small 
% Created: 10/8/99
% Updated: 20/8/99

disp('rb_dPhi called ... rb_dPhi only works for rbm');
[dx,nx]= size(x);
[dc,nc]= size(base.centres);

z= ones(1,nx);

T=zeros(dx,nx,dx+1+nc);

% constant term
%T(:,:,1)= 0;

% linear terms
for i=1:dx,
	T(i,:,i+1)= z;
end;

% radial basis vectors
%z= ones(1,nx);
for i=1:nc
   d=zeros(dx,nx);
   vs= v(base.strategy(i),:);
   is= find(vs>=0);
   vs= 1+vs(is);  
   d(vs,:)= x(vs,:) - base.centres(vs,i)*z;
   if dx>1,
      dd=sqrt(sum(d.*d));
   else,
      dd=sqrt(d.*d);
   end
   ind=find(dd>0);
   d(:,ind)=d(:,ind)./(ones(dx,1)*dd(ind));
%   ind=find(dd==0);
%   d(:,ind)=nan;
   ff= ones(dx,1) * feval( ['d',func{base.func(i)}], dd, base.radii(:,i) ) ;
   T(:,:,dx+1+i) = d .* ff;
end
