function [y_test,y_pred,e_pred,e_driv]= bffreerun(z_test,start,time,noise,bounds,mu)

% [y_test,y_pred,e_pred,e_driv]= bffreerun(z_test,start,time,noise,bounds,mu)
% make a free run prediction of data
% initial point at start 	   (default start=1)
% run for time steps 		   (default time=300)
% add dynamic noise of std noise   (default noise=0.0)
% stop if outside bounds 	   (default bounds=10*std(a))
%
% start is relative to embedding of z_test
% if start<0, then it is relative to z_test in which case
% require -start > window of embedding
%
% y_test = data actually predicted
% y_pred = predictions of y_test
% e_pred = residual prediction errors, so that y_pred= y_test+e_pred
% e_driv = residuals used in simulation
%
% bifurcation parameter included : mu(1) is intial value, and mu(2) is the increment
%if mu=[] or omitted, then defualt values (mu=bfmu)  are used.
%
% M. Small 
% Created: 13/6/99
% Updated: 18/10/99

rb_get_globals
na =nargin;
nout=nargout;

if na<6,
   disp('bifurcation parameter initial value and increment not specified in bffreerun');
   disp('using defualt');
   mu=bfmu;
end;
if na<5
   bounds= inf;
end
if na<4
   noise= 0;
end
if na<3
   time=[];
end
if na<2
   start=1;
end
if na<1,
   z_test=[];
end;


if length(z_test)==0
   z_test=rb_y;
end;

v= rb_embed;
v= v(~isnan(v(:)));
p=-min(v); 
w= max(v) + p;

if p<=0 | w<=0,
   disp('WARNING : nonsense embedding');
end;

v= [min(v)+(0:w)]; 
vx= 1 + v(find(v>=0));

if isempty(time),
   time= length(z_test)-w;
end;

% embed z_test
[xx,yy]= embed(z_test,v+p-1);
[dx,nx]=size(xx);
xx(dx,:)=(0:1:(nx-1))*mu(2)+mu(1);

% work out where to start from
if start<0
   start= -start-w;
end
if (start<0|start+w-1>length(z_test));
   disp('ERROR : start is too close data end');
end;

% z = current state
z= xx(:,start);
zt=xx(:,(start+1):(start+p-1)); 
x_test= xx(vx,start:length(yy));
y_test= yy(start:length(yy));

% zt = history of state
zt= []; 
err=0;
if nout==4, 
  e_driv=[];
end;

mut=mu(1)+mu(2);

for t=1:time
   
   z= rb_iterate(z);
   z(dx)=mut;
   mut=mut+mu(2);
   if sum( z<-bounds | z>bounds) ,disp('exceeded bounds');
      break, 
   end;
   
   if length(noise)==1, 
      err = noise*randn(1);
   elseif length(noise)>1, 
      err = noise(t); 
   end;
   z(1) = z(1) + err;
  
   zt= [ zt z ];
   z=zt(:,t); %formerly dne.
   if nout==4,
      e_driv=[e_driv err];
   end;

end


if min(size(zt))>=1,
   y_pred= zt(1,:);
else
   y_pred=[];
end;

if sum( z<-bounds & z>bounds)
   disp(sprintf('Exceeded bounds at %d',t));
end

t= 1:min(length(y_test),length(y_pred));
y_test=y_test(t);
e_pred= y_pred(t) - y_test(t);








