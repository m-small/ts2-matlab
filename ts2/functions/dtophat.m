function y = dtophat(x,s)

% y = dtophat(x,s)
%
% returns the derivative of the tophat gaussian function with radius (std dev) s(1) and
% curvature s(2).
%
% M. Small 
% Created: 10/8/99
% Updated: 23/9/99


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

%exp(f*abs( )^p) ->0 faster than (abs( )^(p-1))^(-1) so we do the following to avoid NaNs
ex=exp( f*abs(x./s).^p );

ind=find(ex~=0);
y=ex;
if p~=1,
   y(ind) = (f*p/s)*sign(x(ind)).*sign(s).*(abs(x(ind)./s).^(p-1)).* ex(ind);
else,
   y(ind) = (f*p/s).*sign(x(ind)).*sign(s).* ex(ind);
end;

