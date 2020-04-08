function y = dbspline(x,pms)
% Derivative of cardinal B-spline function calculated in bspline.f
%
% function y = dbspline(x,pms)
%
% We use, Nm'(t) = N(m-1)(t) - N(m-1)(t-1).
%
%  INPUT:
%   'x' vector of input arguments
%   'pms' vector of function parameters
%
%  OUTPUT:
%   'y' function values
%
% Notes:
%   'x' is a 1-by-n or n-by-1 vector
%   'y' is 1-by-n
%    y(i) is the evaluation of the function at x(i) with parameters values
%    pms(:,i) (or pms(:) if pms is 2-by-1).
%
% The format of the vector 'pms' is as follows:
%    pms = [a; b; m];
% with
%   'a' dilation parameters (size = [1 n])
%   'b' translation parameters (size = [1 n])
%
% The order of the spline, parameter 'ord', is taken in from a global 
% variable 'SPLINE_ORDER' in the global MATLAB workspace.
%    'ord' order of spline (polynomial of degree ord-1) 
%
% NOTE: If a vector 'x' is given with only one paramter 'ord'. We 
% calculate the output of the B-spline for values 'x' and order 'ord'.
%  If 'pms' (2-by-1) is given as only two parameters, those parameters 
% are taken for all values of 'x'. The functions are currently centred 
% at zero too.
%
% Brian Fleming, 11/6/00.

x = x(:)';
no_params = 2;

if (size(pms,1)*size(pms,2)) ~= 1

   if size(pms,1) ~= no_params
      error('Parameter matrix orientation incorrect!')
   end
   
   if (((size(pms,1)*size(pms,2)) ~= (no_params*length(x)))...
   					   & (size(pms,2) ~= 1))
      error('Incorrect number of parameters to function!')
   end 
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Change global variable
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   global SPLINE_ORDER;
   if length(SPLINE_ORDER) == 0
      error('Spline order not defined!')
   end
   if (SPLINE_ORDER  < 2 | SPLINE_ORDER > 99)
      error('Spline order out of range for derivative!')
   end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   SPLINE_ORDER = SPLINE_ORDER - 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Centre Function

   x = x + pms(1,:)/2;
   
% Note: we have to add this corrective factor 
% because we have changed the SPLINE_ORDER but
% we still require the 'bspline' functions to 
% be centred with the original 'order'.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   y = (1./pms(1,:)).*(bspline(x,pms) - bspline(x-pms(1,:),pms));
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Change global variable
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   SPLINE_ORDER = SPLINE_ORDER + 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
else
   if (pms  < 2 | pms > 99)
      error('Spline order out of range for derivative!')
   end
   pms = pms - 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Centre Function

   x = x + 1/2;
   
% Note: we have to add this corrective factor 
% because we have changed the SPLINE_ORDER but
% we still require the 'bspline' functions to 
% be centred with the original 'order'.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   y = bspline(x,pms) - bspline(x-1,pms);
end





