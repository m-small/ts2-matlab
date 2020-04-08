function [yt,yp,ep,x]=bfpredict(y,mu);

%function [yt,yp,ep]=bfpredict(y,mu);
%
%make predictions on y using global radial basis model - with bifurcation parameter, mu
%ep=yp-yt
%
% to scale mu, take mu=bfmu(mu)
%
% M. Small 
% Created: 7/6/99
% Updated: 13/10/99

rb_get_globals;

if nargin<1,
   y=rb_y;
   mu=bfmu;
   mu=((1:length(y))-1)*mu(2)+mu(1);
end;

v=rb_embed;
v=min(v(~isnan(v))):max(v(~isnan(v)));
[x,yt]=embed(y,v);
[dx,nx]=size(x);
mu=mu(:)';
if length(mu)==1,
   mu=mu*ones(1,nx);
elseif length(mu)<nx,
   disp('ERROR : length of mu and size of X don''t match in bfpredict');
end;
x(dx,:) = mu(1:nx);

yp=rb_image(x)';
ep=yp-yt;


