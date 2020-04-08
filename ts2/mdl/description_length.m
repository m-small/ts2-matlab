function dl = description_length(mss,lambdas,deltas,n,alpha)

% dl=description_length(mss,lambda,delta,dim)
% Description length for a pseudo-linear model
%
% The description length is calculated from the accuracies "delta" of
% the parameters "lambda", using floating point prior.
%
% "mss", the mean sum of squares, is implicitly a
% parameter and its description length is included in the returned value.
%
% mss= mean sum of squares
% lambda= vector of parameters
% delta= vector of parameter accuracies (length= #parms)
% dim= dimension of original data

% Copyright (c) 1994 by Alistair Mees and Kevin Judd.
% Please see the copyright notice included
% in this distribution for full details.
%
%   $Id: description_length.m,v 1.1 1994/11/09 08:54:46 kevin Exp $

if nargin<5
  alpha=0;
end;
%warn against nosense
if any(deltas<0), disp('*** Some deltas are negative ***'); end;

%warning(any(abs(lambdas)<deltas),'*** Some lambdas < deltas ***');

% number of parameters and error
k= length(deltas);

% description length of parameters: !!!THIS ASSUMES DELTAS ARE ABSOLUTE!!!
if k>0
  prc= sum(Lstar(lambdas,deltas));
else
  prc= 0;
end

% dl is prc plus descr length of errors
dl= 0.5*n*(1+log(2*pi*mss)) + prc + alpha*k;
