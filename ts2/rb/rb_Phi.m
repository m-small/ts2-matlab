function X = rb_Phi(x,base,v,func,meth)

% X = rb_Phi(x,base,v,func)
% 
% make a matrix X  at centres base.centre of radius base.radii 
% using embedding strategy base.strategy of v and basis function base.func.
%
% Modified to look out for time embeddings - Infs in v.
% timescale is the affine transformation for time.
%
% M. Small 
% Created: 31/3/98
% Updated: 6/6/99


[dx,nx]= size(x);
[dc,nc]= size(base.centres);

T= [];
z= ones(1,nx);

% constant term
if any('c'==meth),
  T= z';
end;

% linear terms
if any('l'==meth),
  T= [ T x' ];
end;

% radial basis or wavelet basis vectors
if any('r'==meth),
  %radial basis vectors
  %basis functions are \Phi((x-c).^2/r) 
  %for rbm
  [dt,xc]= size(T);
  X= zeros(nx,nc+xc);
  if xc>0,
   X(:,1:xc)= T;
  end
  z= ones(1,nx);
  for i=1:nc
    vs= v(base.strategy(i),:);
    is= find(vs>=0);
    vs= 1+vs(is);
    d= x(vs,:) - base.centres(vs,i)*z;
    if length(vs)==1,
      dd=sqrt(d.*d);
    else,
        dd=sqrt(sum(d.*d));
    end
    X(:,xc+i) = feval( func{base.func(i)}, dd, base.radii(:,i) )';
  end
  %no nonlinear terms.
  if nc<1,
    X= T; 
  end
elseif any('w'==meth),
  %wavelet basis vectors
  %basis functions are \Phi((x-c)/r) 
  %for mdl best wavelet expansion
  [dt,xc]= size(T);
  X= zeros(nx,nc+xc);
  if xc>0,
   X(:,1:xc)= T;
  end
  z= ones(1,nx);
  for i=1:nc
    vs= v(base.strategy(i),:);
    is= find(vs>=0);
    vs= 1+vs(is);
    d= x(vs,:) - base.centres(vs,i)*z;
    %this is the only bit thats different (to rbfs) --- its just
    %quicker to have if outside the for loop.
    if length(vs)==1,
      dd=d;
    else,
      dd=sum(d); %wavelets only REALLY make sense in 1D but...
    end
    X(:,xc+i) = feval( func{base.func(i)}, dd, base.radii(:,i) )';
  end
  %no nonlinear terms.
  if nc<1,
    X= T; 
  end
elseif any('n'==meth),
  %neural net
  %basis functions are \Phi(x.c-r) 
  %for ann
  [dt,xc]= size(T);
  X= zeros(nx,nc+xc);
  if xc>0,
   X(:,1:xc)= T;
  end
  for i=1:nc
    vs= v(base.strategy(i),:); %leave all this embedding stuff in
    is= find(vs>=0);           %but it does no good anyway
    vs= 1+vs(is);              %embedding stuff
    dd= base.centres(vs,i)'*x(vs,:);  %this bit is different
    X(:,xc+i) = feval( func{base.func(i)}, dd-base.radii(1,i) )';
  end
  %no nonlinear terms.
  if nc<1,
    X= T; 
  end
elseif any('d'==meth),
    %use the so-called (n-1)-D doublet basis function
    % the last component of X serves as Kronecker delta - multiplied by
    % the usual basis function evaluated over the others...
    global rb_bandwidth rb_gap
    [dt,xc]= size(T);
    X= zeros(nx,nc+xc);
    if xc>0,
        X(:,1:xc)= T;
    end
    for i=1:nc
    vs= v(base.strategy(i),:); %leave all this embedding stuff in
    is= find(vs>=0);           %but it does no good anyway
    vs= 1+vs(is);              %embedding stuff
    d= x(vs,:) - base.centres(vs,i)*z;
    if length(vs)==1,
      dd=sqrt(d.*d);
    else,
        dd=sqrt(sum(d.*d));
    end
    dd= base.centres(vs,i)'*x(vs,:);  %this bit is different
    X(:,xc+i) = feval( func{base.func(i)}, dd, base.radii(:,i) )'.*( (base.centres(end,i)-x(end,:))<rb_bandwidth/2 | (base.centres(end,i)-x(end,:)+rb_gap*rb_bandwidth)<rb_bandwidth/2)';
  end
  %no nonlinear terms.
  if nc<1,
    X= T; 
  end
    
elseif any('m'==meth),
  %multiplicative 1-D basis functions
  %basis function are \Prod\Phi((x-c)/r)
  %
  [dt,xc]= size(T);
  X= zeros(nx,nc+xc);
  if xc>0,
   X(:,1:xc)= T;
  end
  z= ones(1,nx);
  for i=1:nc
    vs= v(base.strategy(i),:);
    is= find(vs>=0);
    vs= 1+vs(is);
    d= x(vs,:) - base.centres(vs,i)*z;
    Xp=ones(nx,1);
    for j=1:length(vs),
        Xp=Xp.*feval( func{base.func(i)}, d(j,:), base.radii(:,i) )';
    end;
    X(:,xc+i) = Xp;
  end
  %no nonlinear terms.
  if nc<1,
    X= T; 
  end
end;


