function y = wavelet(x,s)

% y = wavelet(x,s)
%
% returns wavelet function of size s
%    X=x/s
%    y=(2X^2-1)*exp^(-X^2)
% this is the wavelet basis function with "relatively low curvature" described by Allingham, et al.
% (the mexican hat function).
%
% M. Small 
% Created: 18/6/99
% Updated: 18/6/99

if nargin<2
  s=1;
end
if length(s)<1
  s=1;
end

s = s(1,:);

x=x./s;

y = (2*x.^2-1).*exp(-x.^2);

