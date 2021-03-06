# ts2-matlab
timeseries modelling via radial basis functions and minimum description length, implemented in MATLAB

This code should be good to go - install all directories somewhere where MATLAB can find them and then place *all* these directories at the top of MATLAB's search path (the code is old and uses some function names that MATLAB has since used - strangely, one of these is trace).

Then, at the command line, suppose you have a time series y,

``` MATLAB
trace('on',1)

Rissanen = 'dl1(mss,a,deltas,dim_x)';

v=genve([-1 0 10 20]) %variable embedding with lag 10, embedding dimension 3, and prediction step 1
func={'gaussian','tophat'}; %basis functions are Gaussians and tophat Gaussians - there are other options
%build a model
buildmodel(y,2000,v,func,Rissanen,3,5,250,3,0,'clr')
% model built on data in vector v, first 2000 points to fit, the rest to test.
% v is the embedding strategies, func the basis functions and penalty the information
% criteria. Repeat 3 times, 5 progressive builds each (this is a greedy algorithm)
% 250 candidate basis functions at each step, go 3 steps past the minimum of DL
% don't apply the local nonlinear optimisation (0->1 otherwise), and include
% Constant, Linear and Radial (i.e. 'clr') terms ('cln' will build a sinlge layer
% network.).
%
% Note - if y is a noise free simulation of a deterministic dynamical system, this
% algorithm could go for some time. Try y+randn(size(y)*std(y)/10, for example
```

Once the model is built it is stored as global variables. It can be applied to
do onestep predictions (```predict```) or long term simulations (```freerun```).
Running ```buildmodel``` will endeavour to improve on the existing model. Sometimes,
this will not be possible because MDL calculations between models are not always
comparable. In such situations, use ```rb_clear_globals``` to remove an existing
model stored in global.

Have fun.

If used for scientific purposes, please cite:
+ M. Small and C.K. Tse. "Minimum description length neural networks for time series prediction." Physical Review E 66 (2002), 066701. pdf
+ M. Small and K. Judd. "Comparison of new nonlinear modeling techniques with applications to infant respiration." Physica D, 117 (1998): 283-298.
+ K. Judd and A. Mees, "On selecting models for nonlinear time series"
Physica D 82, 426 (1995).

Michael Small
upload to github April 2020
