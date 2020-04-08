function [mu]=rb_wobbly(xi,X,v,func,e,meth,alpha);

% function [mu,phi]=rb_wobbly(xi,X,v,func,e,meth,alpha);
%
% evaluation of obj=mu=-phi'*e
%
% see rb_wobble
%
% modified to incorporate rh_wobbly_Phi
%
% M. Small 
% Created: 26/6/99
% Updated: 11/2/02

beta=9;

is= find(v>=0);
vs= 1+v(is);
lvs=length(vs);

base.centres(vs,1)=xi(1:lvs);
base.radii=xi((lvs+1):end);
base.strategy=1;
base.func=1;
phi=rb_Phi(X,base,v,func,meth);

phi=normalize(phi);
mu=-abs(phi'*e); %this is the sensitivity

%the elastic band principle:
if any(meth=='r'),
  %this is reasonable is basis functions are radial
  %else, not.
  [dx,nx]=size(X);
  pen=X(vs,:)-base.centres(vs,1)*ones(1,nx); % distance to data
  if dx>1,                       % square distance to data
    pen=sum(pen.*pen);
  else,
    pen=pen.*pen;
  end;
  sx=std(X(1,:));  
  pen=min(pen)./(sx);            % normalised minimum distance to data
  pen=max(alpha,pen);            % penalise if greater than alpha
  mu=mu./((1-alpha)+pen).^beta;  % penalise for have basis functions away from data
  %actually, there is two parameters: alpha, which specifies when the
  %penalty kicks in; and beta, which specifies how sharp it is.
  
  %disp(num2str((1-alpha)+pen));
  %the elastic band principle:
elseif any(meth=='d'),
%for the twin peak fdoublet, we can do this too
% just don't go in the last direction
  [dx,nx]=size(X);
  vs(vs==dx)=[]; %don't go there
  pen=X(vs,:)-base.centres(vs,1)*ones(1,nx); % distance to data
  if dx>1,                       % square distance to data
    pen=sum(pen.*pen);
  else,
    pen=pen.*pen;
  end;
  sx=std(X(1,:));  
  pen=min(pen)./(sx);            % normalised minimum distance to data
  pen=max(alpha,pen);            % penalise if greater than alpha
  mu=mu./((1-alpha)+pen).^beta;  % penalise for have basis functions away from data
  %actually, there is two parameters: alpha, which specifies when the
  %penalty kicks in; and beta, which specifies how sharp it is.
  
  %disp(num2str((1-alpha)+pen));
end;




