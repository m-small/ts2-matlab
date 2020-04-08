function [ind,x,y]=gettrain(y,v,n,k);

%function [ind,x,y]=gettrain(y,v,n,k);
%
%select the "optimal" training set of size n from the scalar time series y
%when using the embedding v and the k nearest neighbours. ind is the index
%of the selected embedded points such that if [xx,yy]=embed(y,v), then the
%output x=xx(:,ind), y=yy(:) and the model should fit f(x)~y.
%
% Default
% k=10
% n=1000
% v=[-1 0 1 2 3 4]

narg=nargin;

if narg<4,
    k=10;
    if narg<3,
        n=1000;
        if narg<2,
            v=[-1 0 1 2 3 4];
        end;
    end;
end;

%embed
[x,y]=embed(y,v);
[dx,nx]=size(x);

%is n<nx?
n=min(n,nx);

%start
ind=1:n;
i=n+1;

%go!
figure(gcf);
plot(y);
hold on;
p=plot(ind,y(ind),'r.');
set(gcf,'DoubleBuffer','on')
hold off;
while i<nx,    
    %do local constant for now
    rr=rms(x(:,ind)'-ones(n,1)*x(:,i)');
    [dist,near]=sort(rr);
    near=ind(near);
    e1=rms(y(near(1))-y(i));
    
    %guess a point
    ref=ceil(rand(1)*(n-k-1)+k+1);
    ref=near(ref);
    ref=ind(ind==ref);
    
    %a predict/interpolate it too 
    rr=rms(x(:,ind)'-ones(n,1)*x(:,ref)');
    rr(ind==ref)=rms(x(:,i)'-x(:,ref)');
    [dist,near]=min(rr);
    near=ind(near);
    if near==ref,
        near=i;
    end;
    e2=rms(y(near)-y(ref));
    
    %exchange?
    if e2<e1,
        ind(ind==ref)=i;
    end;
    i=i+1;
    
    set(p,'XData',ind);
    set(p,'YData',y(ind));
    drawnow;
end;

%close
set(gcf,'DoubleBuffer','off')

%output
narg=nargout;
if narg>1,
    x=x(:,ind);
    y=y(ind);
end;
