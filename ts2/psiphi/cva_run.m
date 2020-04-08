function [y_test,y_pred,e_pred]= cva_run(z_test,start,time,noise,A,sc,np)

% [y_test,y_pred,e_pred]= cva_run(a,start,time,noise,A,sc,np)
% 
% just like pl_run but for free runs when prediction more than one step into
% the future provided A and sc are given (see cva)
%
% %make a free run prediction of data
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
%
% Modifications:
% Adjusted to work(?) if some of the pl_lags are inf. In this case it assumes
% that they are time coordinates.
% if mu is specifies then it is the (fixed) values of the bifurcation/time 
% parameter.
%
% Modified: in the case of prediciting more than one time step at a time (ie 
% v(:,1)<-1) y_pred and y_test are the approriate predictions and errors.
% Now to predict p steps into the fututre  predictions are made for every
% point. I.e. iterative predictions of the first p points, rather than 
% just iterating one initial point. The result is the composition of 
% independent free run predictions of the first p points under the same 
% psuedo linear model. Hence if the map is not perfect it the model will 
% become eratic rapidly. But what would be better???
%
% Basically, making any free run prediction when predicting more than one step
% into the future is going to give rather pitiful results.
%
% 9/2/96
% 
% If noise='de' it predicts the errors using density estimation of residuals
% from one step predictions of z_test.
% If noise='de',int2str(x) it predicts the errors using (x)% of the residual
% predicted by density estimation (actually (x)% of the deviation of the 
% residual predicted by de from the mean residual predicted by de).
%
% 12/11/96
% 
% added bifurcation parameter mu.
% 9/1/97
% MALS
 
% Copyright (c) 1994 by Kevin Judd.  
% Please see the copyright notice included
% in this distribution for full details.
%
%   $Id: pl_run.m,v 2.1 1995/03/08 09:44:03 kevin Exp $


pl_globals


na =nargin;

if length(z_test)==0
  z_test=pl_y;
end;

v= pl_lags;
p=-v(1); % dne. 
w= v(length(v)) - v(1);

v= [min(pl_embed(:,1))+(0:w)]; 	% used to be v=[-1 0:w-1]; 
        	% this is only OK if only predicting one step ahead - or if not
		% using a variable embedding.
vx= 1 + v(find(v>=0));
  bounds= 10*std(z_test);
  mu=[];

if na<6
  sc=1;
end;
if na<5
  A=1;
end
if na<4
  noise= [];
end
if na<3
  time= length(z_test)-w;
end
if na<2
  start=1;
end
if isstr(noise)
  if strcmp(noise(1:2),'de')
    de=1;
    if length(noise)>2
      noise=eval(noise(3:length(noise)))./100;
    else
      noise=1;
    end;
  else
    disp('WARNING: unknown string specified to pl_run, assuming de.');
    de=1;
  end;
else
  de=0;
end;


% embed z_test
[xx,yy]= embed(z_test,v+p-1);
if any(any(pl_embed==inf)) % a time co-ordinate? - if so add another col
  vv= pl_embed(find((pl_embed>=0 | pl_embed<0) & pl_embed<inf));
  vs= unique(vv);
  [a,b]=size(pl_x); [a,c]=size(xx);
  timescale=(max(pl_y)-min(pl_y))/b;
  if ~isempty(mu),
    xx(a+1,:)=mu;
  else 
    xx(a+1,:)=(1:c).*timescale+min(pl_y);
  end;
end;


% work out where to start from
if start<0
  start= -start-w;
end
fatal(start<0|start+w-1>length(z_test),'start is too close data end');

% z = current state
z= xx(:,start);
zt=xx(:,(start+1):(start+p-1)); % formerly dne.
x_test= xx(vx,start:length(yy));
y_test= yy(start:length(yy));

% zt = history of state
zt= []; %formerly did exist.

for t=1:time
    z= pl_image(z);
    adj=A*(polyn(flipud(z(1:p))',np)./sc)';
    adj=flipud(adj);
    z(1)=adj(1);
  
    zt= [ zt z ];
    z=zt(:,t); %formerly dne.

  %  zt(1,(1:p)+t-1)=adj';
  %  zt(:,(2:p)+t-1)=[];
if 0,
  for i=1:p
    z= pl_image(z);
   
  
    zt= [ zt z ];
    z=zt(:,t+(i-1)); %formerly dne.
  end;
    adj=A*(polyn(zt(1,(1:p)+t-1),np)./sc)';
    zt(1,(1:p)+t-1)=adj';
    zt(:,(2:p)+t-1)=[];
end;
end

if sum( z<-bounds & z>bounds)
  disp(sprintf('Exceeded bounds at %d',t));
end

if min(size(zt))>=1,
 y_pred= zt(1,:);
else
 y_pred=[];
end;

t= 1:min([length(y_test),length(y_pred)]);
e_pred= y_pred(t) - y_test(t);








