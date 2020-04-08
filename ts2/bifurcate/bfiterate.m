function [y_test,y_pred,e_pred,merr,corr]= bfiterate(z_test,mus,time)

% [y_test,y_pred,e_pred,merr,corr]= bfiterate(z_test,mus,time)
% or x_pred=bfiterate(x_pred,mus,time);
%
% make a free run predictions (iterated predictions) of data
% run for time steps 		   (default time=300)
% keep bifurcation parameter value fixed at mus
%
% y_test = data actually predicted - a matrix elementwise, it is comparable to y_pred
% y_pred = predictions of y_test - a matrix, y_pred(time-i+1,:) are the i-step predictions
% e_pred = residual prediction errors, so that y_pred= y_test+e_pred
% merr = mean prediction error - merr(i) is the mean i-step prediction error
% corr = correlation between predicted and actual values - corr(i) is the correlation between 
% actual and i-step predictions.
%
%
% to scale mu, take mu=bfmu(mu)
%
% M. Small 
% Created: 3/7/99
% Updated: 14/10/99

rb_get_globals
na =nargin;
nout=nargout;

if na<2
   time=100;
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
%w= max(v) + p;
w=max(v);

if p<=0 | w<=0,
   disp('WARNING : nonsense embedding');
end;
if p~=1,
   disp('large prediction step - this ain''t going to work');
end;


%v= [min(v)+(0:w)]; 
v=[(0:w)];
vx= 1 + v(find(v>=0));

if isempty(time),
   time= min(100,length(z_test)-w);
end;

% embed z_test
if min(size(z_test))==1,
  xx= embed(z_test,v);
  y_test = embed(z_test((w+2):end),time,1);
else,
   xx=z_test;
   y_test=[];
end;
clear z_test;
[dx,nx]=size(xx);
xx(dx,:)=mus(1:nx);

% work out where to start from

% z = current state
if ~isempty(y_test),
   z= xx(:,1:(nx-time));
   mus=mus(1:(nx-time));
else
   z= xx(:,1:nx);
   mus=mus(1:nx);
end;

% zt = history of state
y_pred = []; 
err = 0;

for t=1:time
   
   z= rb_iterate(z);
   z(dx,:)=mus;
  
   y_pred= [ z(1,:); y_pred];

end

if isempty(y_test);
   y_test=y_pred;
   clear y_pred;
elseif nout>2,
   e_pred= y_pred - y_test;
   if nout>3,
      merr=flipud(rms(e_pred));
      if nout>4,
         corr=corrcoef([y_pred' y_test']);
         corr=corr(1,(time+1):end);
      end;
   end;
end;








