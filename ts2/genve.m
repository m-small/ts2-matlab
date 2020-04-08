function v = genve(m,p,q)

% v=rb_genve(m,p)
% generate variable embedding
% 
% this does it the way I want it to - it also works if m is a vector.
% m if the lags (or maximum lag), and p is the prediction step, q can
% be inf.
% 
% M. Small 
% Created: 21/11/96
% Updated: 6/6/98


if nargin<2
  p= 1;
end;
if nargin<3
  q=0;
end;
if p==inf
  q=inf; p=1;
end;

m=unique(m);
if max(length(m))>1, vect=sort(m(m>=0)); m=sum(m>=0); else, vect=[]; end;

v= [];
k= 1;
for n= 1:m
  v= [ v zeros(k,1); v n*ones(k,1) ];
  k= 2*k;
end;
v(1,:)= []; 				% delete row of zeros
v= v - 1; 				% indices start from 0
id= find(v<0); 				% find non indices
v(id)= NaN*ones(size(id)); 		% convert to Nan
v= sort(v')';
v= [-p*ones(k-1,1) v];

if q==inf
 [av,bv]=size(v); vv=v;
 vv(:,bv+1)=inf*ones(av,1);
 v(:,bv+1)=nan*ones(av,1);
 v=[v; vv;];
end;
 
v= sort(v')';

if length(vect)==m 
   for i=2:length(v(1,:))
     ind=find(~isnan(v(:,i)) & ~isinf(v(:,i)));
     v(ind,i)=vect(v(ind,i)+1)';
   end;
end;


