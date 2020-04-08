function [yt,yp,ep,x]=predict(y);

%function [yt,yp,ep]=predict(y);
%
%make predictions on y using global radial basis model
%ep=yp-yt
%
%
% M. Small 
% Created: 7/6/99
% Updated: 7/6/99

rb_get_globals;

if nargin<1,
   y=rb_y;
end;


v=rb_embed;
v=min(v(~isnan(v))):max(v(~isnan(v)));
[x,yt]=embed(y,v);
yp=rb_image(x)';
ep=yp-yt;


