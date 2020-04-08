function y=rb_image(X);

%function y=rb_image(X);
% 
%calculate the image of columns of X under the global radial basis model. 
%i.e. y=f(X);
%
% M. Small 
% Created: 7/6/99
% Updated: 6/6/99

rb_get_globals;
[d,np]=size(X);

phi=rb_Phi(X,rb_base,rb_embed,rb_functions,rb_method);
%[phi,scale]=normalize(phi); %necessary to avoid ill-conditioning ???
if isempty(rb_basis),
    y=zeros(np,1);
else,
    y=(phi(:,rb_basis)./(ones(np,1)*rb_lambda(:,2)'))*rb_lambda(:,1);
end;
%y=phi(:,rb_basis)*rb_lambda(:,1);

