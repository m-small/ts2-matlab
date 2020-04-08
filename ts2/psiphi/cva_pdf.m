function [yp,img]=cva_pdf(yt,n_obs,n_dyn,n_it,pv,A,Asc,ANp,is,ry)

% [yp,img] = ...
%       cva_pdf(yt,n_obs,n_dyn,n_it,pv,A,Asc,ANp,image_size,[miny maxy]) 
%
% Calculate the pdf of future values for the point yt with dynamic and
% observational noise given by n_dyn and n_obs respectively. Do this using
% n_it nonlinear iterations. pv is the porediction step vector A, Asc, and
% ANp are the linear transformation matrix, the rescaling of that matrix and
% the order of the polynomial interpolation respectively.
%
%


% Copyright (c) 1998 by Michael Small.
%
% Please see the copyright notice included in this distribution
% for full details.
%
%
% File   cva_movie.m
%   $Id$
%
% Created by Michael Small (<watchman@>) on Wed Aug 12 1998
%
% $Log$

na=nargin;

if na<10,
  maxy=[];miny=[];
else
  maxy=ry(2);
  miny=ry(1);
  clear ry
end;
if na<9,
  is=20;
end;
if na<8,
  ANp=1;
end;
if na<7,
  Asc=1;
end;
if na<6,
  A=1;
end;
if na<5,
  pv=1;
end;
if na<4,
  n_it=[];
end; 
if na<3,
  n_dyn=[];
end; 
if na<2,
  n_obs=[];
end; 

if isempty(n_it),
  n_it=50;
end;
if isempty(n_dyn),
  n_dyn=1;
end;
if isempty(n_obs),
  n_obs=1;
end;

[ps,na]=size(A);
[dy,ny]=size(yt);
if ny>1 & dy==1,
  yt=yt';
  [dy,ny]=size(yt);
end;
  
%put noise (uncertaintity on observations) and generate the large
%i.c. matrix.
x=yt*ones(1,n_it);
x=x+randn(dy,n_it)*n_obs;

%iterate these things out;
pl=cva_evals(x,pv,n_dyn);

%and multiply it all together
yp=cva_pred2(pl,A,Asc,ANp);

%yp should be our realisations
if nargout>1,
  %calculate the image if necessary;
  if isempty(miny),
    miny=min(yp(:)); maxy=max(yp(:));
  end;
  ss=(maxy-miny)/(is);
  grd=miny:ss:maxy;
  for i=1:is;
    img(i,:)=sum(yp'<grd(1+i) & yp'>grd(i));
  end;
end;
  

% End of cva_pdf.m





