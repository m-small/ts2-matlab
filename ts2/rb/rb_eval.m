function [y,X] = rb_eval(x)

% [y,X] = rb_eval(x)
% evaluate the global rb model at columns of x
% X is the evluation of the b.f.s at the cols of x.
%
% M. Small 
% Created: 31/3/98
% Updated: 3/7/99


rb_get_globals;

% extract model data
[dc,nc]= size(rb_base.centres); 				% centres tell us dimension

[dx,nx]= size(x);

v= rb_embed;
[nv,dv]= size(v);

X=rb_Phi(x,rb_base,v,rb_functions,rb_method);

ll=length(rb_basis);
[dX,lX]=size(X);
if ll~=length(rb_basis), 
 disp(['WARNING: ll~=lX, ll=',int2str(ll),'; lX=',int2str(lX),'.']);
 %X= X(:,1:ll);
end;

%if length(rb_base(1,:))>1,
X(:,rb_basis)=X(:,rb_basis)./(ones(dX,1)*rb_lambda(:,2)');
%end;
y= X(:,rb_basis) * rb_lambda(:,1);

if nargout>1,
   X=X(:,rb_basis);
end;
