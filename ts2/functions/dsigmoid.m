function y=dsigmoid(x,s);

%function y=dsigmoid(x,s);
%
% returns d/dx of the sigmoid tanh(mx+b)
% where s(1)=m, s(2)=b;
%
%
% M. Small 
% Created: 31/3/98
% Updated: 23/9/99

if nargin<2,
   s=1;
end;
if length(s)<2,
   s(2)=0;
end;

y = (1-tanh( (x+s(2))/s(1) ).^2)/s(1);


