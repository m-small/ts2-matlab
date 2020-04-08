function dl = dl2(mss,lambdas,deltas,n,base,basic,offset,v)

% dl=dl2(mss,lambda,delta,dim,base,basic,offset,v)
% Description length for a pseudo-linear model
%
% The description length is calculated from the accuracies "delta" of
% the parameters "lambda", for each of the basis functions base subject
% to the embedding strategies v, using floating point prior.
%
% penalty for each parameter associated with a basis function is assumed
% to be the same and be the same as the penalty of the coefficient, 
% lambda.
%
% "mss", the mean sum of squares, is implicitly a
% parameter and its description length is included in the returned value.
%
% mss= mean sum of squares
% lambda= vector of parameters
% delta= vector of parameter accuracies (length= #parms)
% dim= dimension of original data
%
% M. Small 
% Created: 31/3/98
% Updated: 12/9/00

if nargin<5
  
end;
%warn against nosense
if any(deltas<0), disp('*** Some deltas are negative ***'); end;

%warning(any(abs(lambdas)<deltas),'*** Some lambdas < deltas ***');

% number of parameters and error
k= length(deltas);

% description length of parameters: !!!THIS ASSUMES DELTAS ARE ABSOLUTE!!!
if k>0
   if any(basic>offset),
    nparam = sum(v(base.strategy,:)'>=0) + ... % #param is #significant coords. in v +
      sum(base.radii~=0);                  % #radii terms for that func
   nparam(basic>offset)=nparam+1;
   end;
   nparam(basic<=offset)=1;
  
  prc= nparam*(Lstar(lambdas,deltas));
else
  prc= 0;
end

% dl is prc plus descr length of errors
dl= 0.5*n*(1+log(2*pi*mss)) + prc;
