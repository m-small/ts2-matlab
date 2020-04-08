function [z_pred,yp,ep,y_data,B,Bsc,C,Csc]=cva_analysis(data,td,ps,emb,bnp,cnp)

% [z_pred,y_pred,e_pred,y_data,B,Bsc,C,Csc] = ...
%                          cva_analysis(data,testdatum,predstep,emb,BNp,CNp);  
%    typically the cva functions may be used like this
% 
% data is the scalar timeseries
% testdatum is separates the modelling and test data (i.e data(1:testdatum)
%  for modelling, the rest for testing)
% predstep is the maximum prediction time-step 
% emb is the predicted values utilise in the nonlinear model.
% BNp and Cnp are the orders of the polynomial interpolation, BNp for B and
%                CNp for C.
%
% is is assumed that a one-step radial basis prediction model resides in the
% global variables (i.e. there exists a pl_model).
%
% Defaults:
% CNp=BNp
% BNp=3
% emb=[1 3 6 9 12 18 24]
% predstep=60
% testdatum=4001

% Copyright (c) 1998 by Michael Small.
%
% Please see the copyright notice included in this distribution
% for full details.
%
%
% File   cva_analysis.m
%   $Id$
%
% Created by Michael Small (<watchman@>) on Wed Aug 19 1998
%
% $Log$

na=nargin;
if na<6
  cnp=[];
end;
if na<5,
  bnp=[];
end;
if na<4,
  emb=[];
end;
if na<3,
  ps=[];
end;
if na<2,
  td=[];
end;



if isempty(bnp),
  bnp=3;
end;
if isempty(cnp),
  cnp=bnp;
end;
if isempty(emb),
  emb=[1 3 6 9 12 18 24]; %don''t knock it --- it works for the Lorenz.
end;
if isempty(ps),
  ps=60; %ditto!
end;
if isempty(td),
  td=floor(ld*0.8);
end;

%sort out the data
ld=length(data);
fit=data(1:(td-1));
test=data(td:ld);

clear data ld

%embed the fit data
disp(']] [z_pred,y_data]=cva_makeit(fit,ps,emb);');
[z_pred,y_data]=cva_makeit(fit,ps,emb);
clear fit

%build the B matrix
disp(']] [B,Bsc,yp,yt,ep,pl]=cva2(z_pred,y_data,bnp,ps);');
[B,Bsc,yp,y_data,ep,z_pred]=cva2(z_pred,y_data,bnp,ps);                

%build the C matrix
e_data=abs(ep);
disp(']] [C,Csc,ep,e_data,eep,pl]=cva2(z_pred,e_data,cnp,ps);');
[C,Csc,ep,e_data,eep,z_pred]=cva2(z_pred,e_data,cnp,ps);
clear e_data eep

%and show it
if trace,
  f1=figure;
  if strcmp(lower(computer),'pcwin'),
    set(f1,'Position',[5 175 800 380]);%windows, windows, windows everywhere.
  else,
    set(f1,'Position',[19 409 1103 433]);
  end;
  set(f1,'Name','Fit data');
  disp(']] cva_movie(yp,ep,y_data,z_pred,min(1000,min(length(y_data),length(test)));');
  cva_movie(yp,ep,y_data,z_pred,min(length(test),length(y_data)));
  close;
end;

%now grab the test data
%and embed it
disp(']] [z_pred,y_data]=cva_makeit(test,ps,emb);');
[z_pred,y_data]=cva_makeit(test,ps,emb);
clear test ps emb

%and make the predictions
disp(']] y_pred=cva_pred2(z_pred,B,Bsc,bnp);');
yp=cva_pred2(z_pred,B,Bsc,bnp);
disp(']] e_pred=cva_pred2(z_pred,C,Csc,cnp);');
ep=cva_pred2(z_pred,C,Csc,cnp);
clear bnp cnp

%and look at the results
if trace
  f1=figure;
  if strcmp(lower(computer),'pcwin'), %windows!!!
    set(f1,'Position',[5 175 800 380]);
  else
    set(f1,'Position',[19 409 1103 433]);
  end;
  set(f1,'Name','Test Data');
  disp(']] cva_slomo(y_pred,e_pred,y_data,z_pred,length(y_data));');
  cva_slomo(yp,ep,y_data,z_pred,length(y_data));
  close
end;

no=nargout;
if no<7,
  clear C Csc
end;
if no<5
  clear B Bsc
end;
if no<4
  clear y_data
end;
if no<3,
  clear ep
end;
if no<2,
  clear yp
end;
if no<1,
  clear z_pred
end;

clear no

% End of cva_analysis.m
