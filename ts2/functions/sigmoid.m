function y=sigmoid(x,s);

%function y=sigmoid(x,s);
%
% returns the sigmoid tanh((x+b)/m)
% where s(1)=m, s(2)=b;
%
%
% M. Small 
% Created: 31/3/98
% Updated: 31/3/98

if nargin<2,
   s=1;
end;
if length(s)<2,
   s(2)=0;
end;

y = tanh( (x+s(2))/s(1) );


