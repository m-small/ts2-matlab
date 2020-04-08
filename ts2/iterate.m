function [y_test,y_pred,e_pred,merr,corr]= iterate(z_test,time)

% [y_test,y_pred,e_pred,merr,corr]= iterate(z_test,time)
% make a free run predictions (iterated predictions) of data
% run for time steps 		   (default time=300)
%
% y_test = data actually predicted - a matrix elementwise, it is comparable to y_pred
% y_pred = predictions of y_test - a matrix, y_pred(time-i+1,:) are the i-step predictions
% e_pred = residual prediction errors, so that y_pred= y_test+e_pred
% merr = mean prediction error - merr(i) is the mean i-step prediction error
% corr = correlation between predicted and actual values - corr(i) is the correlation between 
% actual and i-step predictions.
%
% M. Small 
% Created: 3/7/99
% Updated: 3/7/99

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
xx= embed(z_test,v);

% work out where to start from

% z = current state
z= xx(:,1:(end-time));
y_test = embed(z_test((w+2):end),time,1);

% zt = history of state
y_pred = []; 
err = 0;

for t=1:time
   
   z= rb_iterate(z);
  
   y_pred= [ z(1,:); y_pred];

end

if nout>2,
   e_pred= y_pred - y_test;
   if nout>3,
      merr=flipud(rms(e_pred));
      if nout>4,
         corr=corrcoef([y_pred' y_test']);
         corr=corr(1,(time+1):end);
      end;
   end;
end;








