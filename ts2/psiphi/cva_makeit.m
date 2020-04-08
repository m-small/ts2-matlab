function [y_pred,y_dat,x_dat] =...
    cva_makeit(data,ps,ts,noise)

% [y_pred,y_dat]=cva_makeit(data,ps,ts,noise)  
%
% hack to make the prediction and fit matrices for multimodels.


% Copyright (c) 1998 by Michael Small.
%
% Please see the copyright notice included in this distribution
% for full details.
%
%
% File   cva_makeit.m
%   $Id$
%
% Created by Michael Small (<watchman@>) on Mon Aug 17 1998
%
% $Log$

if nargin<4,
  noise=0;
end;

global pl_embed
v=pl_embed(~isnan(pl_embed));
maxembed=max(v);
minembed=-min(v);
ps=max(ps,minembed);

[x,y]=embed(data,-ps:1:maxembed);                  
[dx,nx]=size(x);
x_dat=x;

%multiple model predictions go here
if 1,
  mts=max(ts);
  global pl_pred_vect
  pl_pred_vect=-ts+1;
  timer('s');
  y_pred=[];y1=x(1,:)';
  for i=0:mts,
    timer(i/mts);
    if any(i==ts),
      y_pred=[y_pred y1];
    end;
    y1=pl_eval(x);
    x(2:dx,:)=x(1:(dx-1),:);
    x(1,:)=y1'+randn(1,nx)*noise;
  end;
  timer('f');
elseif 0,
  pl_load lorenzes_3  
  y3=pl_eval(x);           
  pl_load lorenzes_6 
  y6=pl_eval(x);
  pl_load lorenzes_9 
  y9=pl_eval(x);
  pl_load lorenzes_12
  y12=pl_eval(x);
  y_pred=[y3 y6 y9 y12];
elseif 0,
  lx=floor(length(x)/3)*3;
  y1=pl_eval(x(:,1:3:lx));
  y2=pl_eval(x(:,2:3:lx));
  y3=pl_eval(x(:,3:3:lx));
  y_pred=[y1 y2 y3];
  y=y(:,1:3:lx);
end;
y_dat=y;            
y_dat=flipud(y_dat);

% End of cva_makeit.m



