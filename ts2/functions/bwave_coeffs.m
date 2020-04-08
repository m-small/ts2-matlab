function q = bwave_coeffs(m)
% Calculates B-wavelet coefficents
%
% function q = bwave_coeffs(m)
% 
% We calculate the B-wavelet coefficients according to:
%  "An Introduction to Wavelets", pp.184, C. Chui.
%
%  INPUT:
%   'm' order of B-spline
%
%  OUTPUT:
%   'q' coefficients
%
% Brian Fleming, 12/6/00.

if m < 1 | m >= 50, error('Invalid Spline Order!'), end

q = [];

for n = 0:3*m-2
   co_sum = 0;
   for l = 0:m
     co_sum = co_sum+(factorial(m)/(factorial(l)*factorial(m-l)))*...
                                                  normbspline(n+1-l,2*m);
   end
   q = [q co_sum*((-1)^n)/(2^(m-1))];
end

% filename = 'bwave_m';
% save([filename int2str(m)],'q')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y = normbspline(x,m)
% bspline function, not centred

y = 0;

% characteristic function

if m == 1
   if x>= 0.0 & x < 1.0
      y = 1.0;
   else
      y = 0.0;
   end
end

% higher order

if m >= 2 & m < 100

   for k = 1:m-1
      a(k) = 0.0;
      x1 = x - k + 1;
      if x1 >= 0.0 & x1 < 1.0
         a(k) = x1;
      end
      if x1 >= 1.0 & x1 < 2.0
         a(k) = 2 - x1;
      end
   end
   
   for p = 1:m-2
      for q = 1:m-1-p
         a(q) = ((x-q+1)*a(q)+(p+q+1-x)*a(q+1)) / (p+1);
      end
   end
   
   y = a(1);
   
end
      


