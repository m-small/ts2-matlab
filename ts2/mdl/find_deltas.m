function deltas = find_deltas(S,deltas,trace)

% deltas=find_deltas(S,deltas)
% Find soln of S*deltas=1./deltas

% Copyright (c) 1994 by Kevin Judd.  
% Please see the copyright notice included
% in this distribution for full details.
%
%   $Id$

%rb_globals

if nargin<3
  trace= 0;
end
if nargin<2
  deltas= [];
end
if length(deltas)>max(size(S)),
  deltas=deltas(1:max(size(S)));
end;

max_dd = 1e-4;
max_err = 1e-4;
min_delta = 1e-6;
max_count=1e2;
count=0;

% First try a subspace restricted problem
% u & w define the subspace, x the linear combination, d are deltas
if isempty(deltas)
   
  [m,n]= size(S);
  u= ones(n,1);
  w= S\u;
  q= S*u;
  uq= [ u q ];
  wu= [ w u ];
  x= [0;1];
  d= wu*x;
  err= uq*x - 1 ./d;
  if trace,disp(sprintf('1st Guess dl= %g',d'*S*d-sum(log(d)))); end;
  dd= 2*max_dd;
  
  % Use Newton's method to solve uq*x- 1 ./ wu*x = 0 by least squares
  while any(abs(dd)>max_dd) & count<max_count
    count=count+1;
    d2= d.*d;
    dd= - [ u+w./d2  q+u./d2 ]\err; 	 % Newton step
    x= x + dd;
    d= wu*x;
    bad= find(d<min_delta); 		       % avoid negative precisions
    d(bad)= min_delta*ones(size(bad));
    err= uq*x - 1 ./d;
  end
  deltas= d;
end

% The following just does Newton but projects to avoid negatives:
% We solve S*deltas-1./deltas==0 using deriv S+1./deltas.^2
deltas = abs(deltas);
err = S*deltas - 1 ./ deltas;
dd = 2*max_dd;

if trace,disp(sprintf('good guess dl= %g',deltas'*S*deltas-sum(log(deltas)))); end;

count=0;
best_deltas=deltas;
%best_err=RMS(err);
%last_deltas=best_deltas; last_err=best_err;
%figure
while any(abs(dd)>max_dd) | any(err>max_err);

  % Newton step:
  SS=(S+diag(1./deltas.^2));
    
  if rcond(SS)<eps;
     disp(['S+diag(1./deltas.^2) close to singular, rcond=',num2str(rcond(SS)),'. Normalizing.']);
     [SS,scale]=normalize(SS);
     dd= -SS\err;
     dd=dd./scale';
  else;
     dd = -(S+diag(1./deltas.^2))\err;
  end;
  deltas = deltas + dd;
  
  % Project any negative deltas back into the real world:
  bad = find(deltas<min_delta);
  deltas(bad) = min_delta*ones(size(bad));
  
  err = S*deltas - 1 ./ deltas;

  count=count+1;
  if count>max_count, 
     disp('WARNING : find_deltas taking too long'); 
     deltas=best_deltas; 
     break; 
  end;
end

if trace,disp(sprintf('final guess dl= %g',deltas'*S*deltas-sum(log(deltas)))); end;

