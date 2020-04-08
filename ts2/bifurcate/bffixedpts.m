function [X0,lam,v,range]=bffixedpts(y,mu);

%function [x0,lam,v]=fixedpts(y,mu);
%
%Compute a (fixed point) bifurcation diagram.
%mu is assumes to be unscaled
%
% M. Small 
% Created: 20/8/99
% Updated: 20/8/99

if nargin<2,
  mu=[];
end;
if nargin<1,
  y=[];
end;

rb_get_globals;
if isempty(y),
  y=rb_y;
end;
[dy,ly]=size(y);
if ly==1,
   y=y';
   [dy,ly]=size(y);
end;
if isempty(mu),
   mu=1:ly;
end;

mu=bfmu(mu);

extra=0.1;
res=100;

%diag of bounding box of data
minx=min(rb_x(:));
maxx=max(rb_x(:));
rx=maxx-minx;
maxx=maxx+extra*rx;
minx=minx-extra*rx;
rx=rx*(1+2*extra);
ss=rx/res;
range=minx:ss:maxx;
[dx,nx]=size(rb_x);
a=ones(dx,1);
x=a*(range);

nout=nargout;

wait=waitbar(0,'bffixedpts : working...');
rmu=mu(end)-mu(1);
k=1;
for m=mu,
  waitbar((m-mu(1))/rmu,wait);
  
  %eval f
  x(dx,:)= m*ones(1,res+1);
  f=rb_eval(x)';



  %plot
  if trace,
    cla;
    plot(range,f-range);title('fixedpts');
    grid;drawnow;hold on;
  end;

  %zeros
  i0=find(diff(sign(f-range)));

  %refine it with bisection method --- NOT ! do something smarter.
  tol=1e-14;
  tol=tol*ss;
  x0=[];
  for i=i0,
    lowx=range(i);lowf=f(i)-lowx;
    highx=range(i+1);highf=f(i+1)-highx;
    while (highx-lowx)>tol,
      midx=(highx+lowx)/2;
      ax=a*midx;
      ax(dx,:)=m;
      midf=rb_eval(ax)-midx;
      if sign(midf)==sign(highf),
         highx=midx;
         highf=midf;
      elseif sign(midf)==sign(lowf),
         lowx=midx;
         lowf=midf;
      else,
         lowx=highx;
      end;
   end;
   x0=[x0 midx];
   if trace,

      disp(['Fixed point at ',num2str(midx)]);
      plot(midx,0,'ro');drawnow;
   end;
   
  end;

  %now calculate eigenvect/val.
  n=length(i0);
  fill=[eye(dx-1) zeros(dx-1,1)];
  if nout>1,
    df=rb_deval(a*x0);
    for i=1:n,
      DF=[df(i,:);fill];
      [u,d]=eig(DF);
      lam(:,i)=diag(d);
      v(:,:,i)=u;
    end;
  end;
  
  X0{k}=x0;
  k=k+1;
end;

close(wait);


