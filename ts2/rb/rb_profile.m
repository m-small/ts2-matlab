function [pr,cdf]=rb_profile(X,np);

%function [pr,cdf]=rb_profile(X,np);
%
%provide a profile of the evluation of the rb model basis functions accross the attractor (rb_x) 
%or X.
%
%cdf(p,i) is the proportion of points on the attractor with an evaluation of basis function i 
%in the range [pr(0,i),pr(p,i)].
%
%plot(pr) plots the evaluation of the basis functions accross the attractor. plot(cdf,pr) does 
%the same thing, normalized by the density of points on the attractor.
%
%evaluate at np points (default np=100).
%
% M. Small 
% Created: 5/8/99
% Updated: 5/8/99

rb_get_globals;

if nargin<2,
   np=[];
end;
if nargin<1,
   X=[];
end;

if isempty(X),
   X=rb_x;
end;
if isempty(np),
   np=100;
end;

[dc,nc]=size(rb_base.centres);
[dX,nX]=size(X);

if dX~=dc,
   disp('WARNING : Sizes don''t match in rb_profile');
end;

np=np-1;
z= ones(1,nX);
for i=1:nc
  vs= rb_embed(rb_base.strategy(i),:);
  is= find(vs>=0);
  vs= 1+vs(is);
  d= X(vs,:) - rb_base.centres(vs,i)*z;
  if length(vs)==1,
    dd=sqrt(d.*d);
  else,
	 dd=sqrt(sum(d.*d));
 end
 mind=min(dd);maxd=max(dd);
 step=(maxd-mind)/np;
 cdf(:,i)=cumsum(hist(dd,mind:step:maxd))'./nX;
 pr(:,i) = feval( rb_functions{rb_base.func(i)}, mind:step:maxd, rb_base.radii(:,i) )';

  
end
  

   