function [y,X] = rb_deval(x)

% [y,X] = rb_deval(x)
% evaluate the derivative of the global rb model at columns of x
%
% y is the derivative of the global model w.r.t each component of x.
% X is the evluation of the b.f.s at the cols of x.
%
% This and rb_dPhi ate coded to evaluate many data points at one time,
% but for a large number of data points this does hog A LOT of memory.
%
% M. Small 
% Created: 10/8/99
% Updated: 11/8/99


rb_get_globals;

% extract model data
[dc,nc]= size(rb_base.centres); 				% centres tell us dimension

[dx,nx]= size(x);

v= rb_embed;
[nv,dv]= size(v);


X=rb_dPhi(x,rb_base,v,rb_functions);

for i=1:length(rb_basis),
   X(:,:,rb_basis(i)) = X(:,:,rb_basis(i))/rb_lambda(i,2);
end;

if nx>1 ,%& nc>1,
   for j=1:dx,
      y(:,j) = squeeze(X(j,:,rb_basis)) * rb_lambda(:,1);
   end;
else
   y = squeeze(X(:,:,rb_basis)) * rb_lambda(:,1);
   y = y';
end;

if nargout>1,
   for i=1:length(rb_base),
      X(:,:,rb_basis(i)) = X(:,:,rb_basis(i))*rb_lambda(i,2);  
   end;   
   X = X(:,:,rb_basis);
end;

