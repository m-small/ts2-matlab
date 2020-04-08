function [x0,lam,v,range]=fixedpts;

%function [x0,lam,v]=fixedpts;
%
%find fixed points x0, eigenvalues lam and eigenvectors v of a global radial basis model.
%
% M. Small 
% Created: 20/8/99
% Updated: 20/8/99

nout=nargout;

extra=0.1;
res=100;

rb_get_globals;

minx=min(rb_x(:));
maxx=max(rb_x(:));
rx=maxx-minx;
maxx=maxx+extra*rx;
minx=minx-extra*rx;
rx=rx*(1+2*extra);
ss=rx/res;
range=minx:ss:maxx;


[dx,nx]=size(rb_x);

%diag of bounding box of data
a=ones(dx,1);
x=a*(range);
%eval f
f=rb_eval(x)';

%plot
if trace,
   plot(range,f-range);title('fixedpts');
   grid;drawnow;hold on;
end;

%zeros
i0=find(diff(sign(f-range)));

%refine it with bisection method
tol=1e-7;
tol=tol*ss;
x0=[];
for i=i0,
   lowx=range(i);lowf=f(i)-lowx;
   highx=range(i+1);highf=f(i+1)-highx;
   while abs(highx-lowx)>tol,
      midx=(highx+lowx)/2;
      midf=rb_eval(a*midx)-midx;
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
if n>0,
  
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
else,
  lam=[];
  v=[];
end;

range=f-range;

if trace,
   disp('Done.');
   hold off;
end;


