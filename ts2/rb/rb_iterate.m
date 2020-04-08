function X=rb_iterate(X,n);

%function Xn=rb_iterate(X,n);
% 
%calculate the n-th iterate of columns of X under the global radial basis model. 
%i.e. Xn=f^(n)(X);
%
% M. Small 
% Created: 7/6/99
% Updated: 6/6/99

if nargin<2,
   n=1;
end;

rb_get_globals;

[d,np]=size(X);

if isempty(rb_basis),
    X=zeros(d,np);
else,

    for i=1:n,
   phi=rb_Phi(X,rb_base,rb_embed,rb_functions,rb_method);
 %  [phi,scale]=normalize(phi); %necessary to avoid ill-conditioning  --- YES!!

    y=(phi(:,rb_basis)./(ones(np,1)*rb_lambda(:,2)'))*rb_lambda(:,1);
%   keyboard;
   X=[y';X(1:(d-1),:);];
end;
end;

