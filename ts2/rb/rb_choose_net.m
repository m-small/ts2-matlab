function [base,phi]= rb_choose_net(x,y,func,v,nc)

% [base,X]= rb_choose_net(x,y,r,e,f,s,v)
% generate centres  for a neural net model
% NB: for a net, the centres are themselves weights,
% it seems reasonable that these weights should be in k[-1/a,1/a]
% where a=min(abs(max(x)),abs(min(x)))
%
% defaults
%   func= 'gaussian';
%   nc=300;
%
% M. Small 
% Created: 3/3/02
% Updated: 4/3/02

[d,nx]= size(x);

fact=1;
datalim=fact*min([abs(max(x'));abs(min(x'))]);
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


% generate centre
base.centres=(2*rand(d,nc)-1).*(datalim'*ones(1,nc));

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

   %scale radii so that they are comparable to centres
   dataspread=dataspread*std(base.centres(:,i)'*x); 
   datacent=mean(base.centres(:,i)'*x); 
   for j=1:nr,
      if isnan(rng(2,j)),
         if isnan(rng(1,j)),
            r=[r; randn(1,nfi)*dataspread+datacent];
         else
            r=[r; randn(1,nfi)*rng(1,j)*dataspread+datacent];
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




