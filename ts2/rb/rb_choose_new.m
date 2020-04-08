function [base,phi]= rb_choose_new(x,y,func,v,nc)

% [base,X]= rb_choose_new(x,y,r,e,f,s,v)
% generate centres using distribution of y
%   y = [] = uniform distribution
% 
% defaults
%   func= 'gaussian';
%   nc=300;
%
% M. Small 
% Created: 7/3/99
% Updated: 7/3/99

[d,nx]= size(x);

dataspread=std(x(1,:));

if nargin<5,
   nc=300;
end;
if nargin<4
  v= 0:d-1;
end
if nargin<3 
  func= 'gaussian';
end
if isempty(func)
   func='gaussian';
end;

ns=length(v(:,1));
if iscell(func),
   nf=length(func);
else,
   nf=1;
   f{1}=func;
   func=f;clear f
end;

% sort errors
if isempty(y)
  w= (0:nx-1)/nx;
else
  w= abs(y);
  w= (w-min(w))/(max(w)-min(w));
end;
[w,wi]= sort(-abs(w));
w= -w;
wm=max(w);

% generate centre
c= zeros(d,nc);
for i=1:nc
  rr=rand*wm;
  j= max(find(w>rr));
  c(:,i)= x(:,wi(j));
end;
base.centres=c;

%generate embedding strategies
base.strategy = floor(rand(nc,1)*ns)+1;

%choose basis functions
base.func= floor(rand(nc,1)*nf)+1;

%choose radii parameters
mp=0;
for i=1:nf,
   mp=max(mp,rb_param(func{i}));
end;
mp=max(mp,1);
base.radii=zeros(mp,nc);
for i=1:nf,
   [nr,rng]=rb_param(func{i});
   r=[];
   ni=find(base.func==i);
   nfi=length(ni);
      
   for j=1:nr,
      if isnan(rng(2,j)),
         if isnan(rng(1,j)),
            r=[r; randn(1,nfi)*dataspread];
         else
            r=[r; randn(1,nfi)*rng(1,j)*dataspread];
         end;
      else,
         r=[r; rand(1,nfi)*(abs(rng(2,j)-rng(1,j)))+min(rng(:,j));];
      end;
   end;
   if nr>0,
      base.radii(1:nr,ni)=r;
   end;
end;

   
% construct X matrix
if nargout>1
  phi= rb_Phi(x,base,v,func);
end




