% script demonstrate  
%
% script to demonstrate some of the "cva"/"future prediction with
% error-bars" modelling results
%
% demonstrate with 
%        y=B[x z] + e
%        x=F(z)
% instead of
%        y=Bx + e
%        x=F(z)

% Copyright (c) 1998 by Michael Small.
%
% Please see the copyright notice included in this distribution
% for full details.
%
%
% File   demonstrate.m
%   $Id$
%
% Created by Michael Small (<watchman@>) on Tue Sep  1 1998
%
% $Log$

%setup path
%directories cva, cva/other and cva/models must be on matlab path 
%(ahead of matlab's functions)
% E.g.
% path('/u/watchman/matlab/cva/:/u/watchman/matlab/cva/other:/u/watchman/matlab/cva/models',path); 

%
%args to cva_analysis are as follows
% [z_pred,y_pred,e_pred,y_data,B,Bsc,C,Csc] = ...
%                     cva_analysis(data,testdatum,predstep,emb,BNp,CNp);  
% data is the scalar timeseries
% testdatum is separates the modelling and test data (i.e data(1:testdatum)
%  for modelling, the rest for testing)
% predstep is the maximum prediction time-step 
% emb is the predicted values utilise in the nonlinear model.
% BNp and Cnp are the orders of the polynomial interpolation, BNp for B and
%     CNp for C. (B is for the prediction and C is for the error bars).
%

%demonstration with a tiny noise level Lorenz
%(Gaussian additive noise with std dev =0.1.
%nonlinear model captures the dynamics.
load smallstep_workspace
pl_load smallstep_model
[z_pred,y_pred,e_pred,y_data,B,Bsc,C,Csc] = ...
    cva_analysis3(y,4001,90,[1 2 3 6 9 12 18 24 36],3,3);
testdatum=4001; predtime=90; predvect=[1 2 3 6 9 12 18 24 36]; BNp=3; CNp=3;

%demonstration with a moderate noise level Lorenz
%(Gaussian additive noise with std dev = 1.
%nonlinear model exhibts a periodic orbit.
load onestep_workspace
pl_load onestep_model
[z_pred,y_pred,e_pred,y_data,B,Bsc,C,Csc] = ...
    cva_analysis3(y,4001,60,[1 2 3 6 9 12 18 24 36],3,3);
testdatum=4001; predtime=60; predvect=[1 2 3 6 9 12 18 24 36]; BNp=3; CNp=3;

%demonstration with a large noise level Lorenz
%(Gaussian additive noise with std dev = 2.
%nonlinear model exhibits a periodic orbit.
load bigstep_workspace
pl_load bigstep_model
[z_pred,y_pred,e_pred,y_data,B,Bsc,C,Csc] = ...
    cva_analysis3(y,4001,60,[1 2 3 6 9 12 18 24 36],3,3);
testdatum=4001; predtime=60; predvect=[1 2 3 6 9 12 18 24 36]; BNp=3; CNp=3;

%dynamic and observation noise Lorenz
%(Both Gaussian, dynamic noise on all three components of the orignal
%system, dyn. & obs. noise std 0.5, and 1 respectively).
%nonlinear model is periodic orbit
load dyn1step_workspace
pl_load dyn1step_model
[z_pred,y_pred,e_pred,y_data,B,Bsc,C,Csc] = ... 
    cva_analysis3(y,4001,60,[1 2 3 6 9 12 18 24 36],3,3);
testdatum=4001; predtime=60; predvect=[1 2 3 6 9 12 18 24 36]; BNp=3; CNp=3;


return

%Finally, and perhaps not so impressive:
%moderate obs noise (As above), and MDL applied to the matrix coefficient
%selection. (this may take a while to execute...).
load onestep_workspace
pl_load onestep_model
[z_pred,y_pred,e_pred,y_data,B,Bsc,C,Csc] = ...
    cva_analysis2(y,4001,60,[1 2 3 6 9 12 18 24 36],3,3);




% End of demonstrate.m




