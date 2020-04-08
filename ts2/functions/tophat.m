function y = tophat(x,s)

% y = tophat(x,s)
%
% returns tophat gaussian function with radius (std dev) s(1) and
% curvature s(2).
%
% M. Small 
% Created: 31/3/98
% Updated: 31/3/98


if nargin<2
   s=1;
end;
if length(s)<2
   p=2;   
else
   p=abs(s(2));
   s=s(1);
end


f=(1-p)/p;

y = exp( f * abs(x ./ s).^p );


