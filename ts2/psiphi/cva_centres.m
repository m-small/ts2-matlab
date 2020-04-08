function [y_pred,y_dat] =...
    cva_centres(data,ps,dontneetit)

% [y_pred,y_dat]=cva_centres(data,ps)  
%
% hack to make the prediction and fit matrices for multimodels.
% the predictions are the centres of the radial basis model stored in the
% global variables.

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

if nargin>2,
  disp('Third argument to cva_centres not used');
end;

pl_globals;
v=pl_embed(~isnan(pl_embed));
maxembed=max(v);
minembed=-min(v);
ps=max(ps,minembed);

[x,y]=embed(data,-ps:1:maxembed);                  
[dx,nx]=size(x);


%multiple model predictions go here
if pl_strategy==1
  X = pl_X( x, pl_method, pl_centres, pl_radii, pl_func );
elseif isempty(pl_ellipse),
  X = pl_X( x, pl_method, pl_centres, pl_radii, pl_func, pl_strategy, v);
else
  X = pl_X( x, pl_method, pl_centres, pl_radii, pl_func, pl_strategy, v , pl_ellipse);
end;
pll=pl_lambda;pll(1)=0;

global pl_pred_vect
pl_pred_vect=-(1:sum(ppl~=0));

X=X(:,pll~=0);
y_pred=X;
y_dat=y;            
y_dat=flipud(y_dat);


% End of cva_makeit.m



