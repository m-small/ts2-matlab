function y = gaussian(x,s)

% y = gaussian(x,s)
%
% returns gaussian function with std dev s
%
% M. Small 
% Created: 31/3/98
% Updated: 31/3/98

if nargin<2
  s=1;
end
if length(s)<1
  s=1;
end

s = s(1,:);

f = -0.5;

y = exp( f * abs(x ./ s).^2 );


