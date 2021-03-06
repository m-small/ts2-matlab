function y = morlet(x,s)

% y = morlet(x,s)
%
% returns morlet wavelet function of size s
%    X=x/s
%    y = cos(2*pi*x).*exp(-(2*(pi/z0)^2)*x.^2) - exp(-z0^2/2-(2*(pi/z0)^2)*x.^2);
%
% defualt z0=5 (or 2nd arguement is [s; z0;])
%
% M. Small 
% Created: 18/6/99
% Updated: 18/6/99

if nargin<2
   s=1;
   z0=5;
end
[ds,ns]=size(s);
if max(ds,ns)<1
   z0=5;
   s=1;
elseif ds>1,
   z0=s(2,:);
   s=s(1,:);
else
   z0=5;
   s=s(:)';
end

x=x./s;

y = cos(2*pi*x).*exp(-(2*(pi/z0)^2)*x.^2) - exp(-z0^2/2-(2*(pi/z0)^2)*x.^2);

