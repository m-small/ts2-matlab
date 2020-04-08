function L = Lstar(x,p)

% L = Lstar(n) and Lstar(x,p)
% Universal prior encoding of positive integers and floats
% One arguement is an integer, Two args a float and precision
%
% Lstar(n)= log(2.865)+(log + loglog + logloglog + ...)(n)

% Copyright (c) 1994 by Alistair Mees and Kevin Judd.
% Please see the copyright notice included
% in this distribution for full details.
%
%   $Id: description_length.m,v 1.1 1994/11/09 08:54:46 kevin Exp $

if any(isinf(x)),
  %it's hopeless, just give up.
  L=inf; 
  disp('WARNING: infinite description length in Lstar - we''ve got an irrelevant centre');
end;

x= abs(x);

if nargin==1
  L= log(2.865);
  x= 1 + ceil(x); 			    % must code for zero too
  while any(x>1);
    x= log(max(1,x));
    L= L + x;
  end
else
%just to make sure it doesn't EVER crash.
  if any(isnan(p)),
    disp('WARNING: NaN detected in Lstar');
    p(isnan(p))=eps*ones(size(find(isnan(p))));
  end;
  p= abs(p);
  % fractional part and sign
  L= Lstar(2*ceil(x./p));
  % exponent and its sign
  L= L + Lstar(ceil(abs(2*log(max([x;1./x])))));
end


