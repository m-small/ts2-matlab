function y = cva_image(xdata,A,sc,np)

% y = cva_image(x,A,sc,np)
% n-th image of columns of x under global pl model
%
% The dimension of x must equal either the window + parms + pred - 1
% 
% default n= 1

% Copyright (c) 1994 by Kevin Judd.  
% Please see the copyright notice included
% in this distribution for full details.
%
%   $Id: pl_image.m,v 2.1 1995/03/08 10:17:26 kevin Exp $
%
% Modification:
% Had to make it so it adds one to the window if there is a time coord.
% See pl_run.
% Now time coord bit actually seems to work.
%
% 1/12/95.


pl_globals;

na=nargin;

  n= 1;
if na<4
  cp=1;
end;
if na<3
  sc=1;
end;
if na<2
  A=1;
end


[dim,k]=size(xdata);
if dim==1 & k>1
  xdata=xdata';
  dim= k;
end


v= pl_lags;
d= length(v);
p= -v(1); 				% prediction
w= v(d); 				% window
v= p + v(2:d);% this is what it used to be but it doesn't work with a var. emb.
%predicting more than one step into the future.
%v= sort(unique(pl_embed)); v=v(~isnan(v)); v=v(~isinf(v));				% indices of required embedded data
d= d-1; 				% dimension
q= pl_parms; 				% parameters

%Do we want a time co-ord?
if any(any(find(pl_embed==inf))) % a time co-ordinate? 
  [a,b]=size(pl_x);              % if so we need a timestep.
  timescale=(max(pl_y)-min(pl_y))/b;
  time=w+p+1; % might as well have a place marker for this as well.
else,
  timescale=0;
end;


% keep these coords for new image
if n==1 & (dim==d | dim==d+q) & all(diff(diff(pl_lags))==0)
  v= 1:d; 				% uniform embedding
  pred= p:dim;
  keep= 1:d-1; 				% keep these coords for new image
  parm= d+(1:q); 			% indices of parameters
else
  pred= p:dim;
 %keep= 1:w;% 1:w+p-1; 			% keep these coords for new image
  keep= 1:(w+p-1);%
  parm= p+w+(1:q); 			% indices of parameters
end;
  
% does xdata contains parameter values?
if dim == d+q 				
  v= [ v parm]; 			% indices of all required components
  keep= [keep parm]; 			% retain for next image
end;

y= [];

for x0=xdata
  x= x0;
  for i=1:n;
    xk=[];
    for j=1:p
      xk=[xk; x(keep+(j-1));];
      z(j,:) = pl_eval( x(pred) );
    end;
    keyboard;
    z=A*(polyn(z',np)./(ones(length(z(1,:)),1)*sc))';z'
    if timescale
    	x= [ z; x(keep); x(time)+timescale];
    else
        x=[ z; xk];
    end;
  end
  y= [ y x ];
end;


