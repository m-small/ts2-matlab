function y_pred =...
    cva_evals(x,ts,noise)

% y_pred=cva_evals(xt,ts,noise)  
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

[dx,nx]=size(x);
    
%global pl_embed
%v=pl_embed(~isnan(pl_embed));
%maxembed=max(v);
%minembed=-min(v);
%ps=max(ps,minembed);

%[x,y]=embed(data,-ps:1:maxembed);                  
%[dx,nx]=size(x);

%multiple model predictions go here
  mts=max(ts);
%  global pl_pred_vect
%  pl_pred_vect=-ts+1;
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

%y_dat=y;            
%y_dat=flipud(y_dat);


% End of cva_makeit.m



