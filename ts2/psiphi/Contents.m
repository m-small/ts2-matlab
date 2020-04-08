% CVA inspired two layer iterated radial basis/polynomial modelling
%
% We've gotta come up with a name for this stuff
%    y = AG(z)+err
%    z = [F(x) F(F(x)) ... F(...F(x)...)]
% F - radial basis model, G set of basis functions (polynomials or whatever)
%
%cva.m                ---  build 2 layer model (beta release)
%cva2.m		      ---  build 2 layer model (use this version)
%cva_analysis.m	      ---  do all the necessary calculation (standard)
%cva_analysis2.m      ---  do all the necessary calculation (apply MDL to A)  
%cva_analysis3.m      ---  do all the necessary calculation (don't bother)
%cva_analysisf.m      ---  do all the necessary calculation (Fourier
%                                             decomposition - don't bother)
%cva_centres.m	      ---  predictions from radial basis funcs
%cva_evals.m	      ---  predictions (z) from 1-step pl_timeseries model
%cva_image.m	      ---  n-step predictions from pl_timeseries model
%cva_iter.m	      ---  hack for n-step predictions from pl_timeseries model
%cva_makeit.bad.m     ---  
%cva_makeit.m	      ---  embed and do pl_predictions for cva type model
%cva_movie.m	      ---  draw some pretty pictures to show long term preds
%cva_movie.old.m      ---  
%cva_pdf.m	      ---  calculate pdf (for multiple predictions)
%cva_pred.m	      ---  make predictions for cva/2-layer model 
%cva_pred2.m	      ---  make predictions for cva/2-layer model (see 
%                                           cva_analysis2) 
%cva_pred4.m	      ---  make predictions for cva/2-layer model (for 
%                                           radial basis not polynomial funcs)
%cva_predf.m	      ---  make predictions for cva/2-layer model (see
%                                           cva_analysisf) 
%cva_reddot.m	      ---  calculate pdf of futures (red dots) for given noise
%cva_run.m	      ---  iterated 1-step predictions of cva/2-layer model
%cva_sensitive.m      ---  like cva.m, but apply full mdl (see cva_analysis2)
%cva_slomo.m	      ---  interactive cva_movie
%cva_view.m	      ---  cva_slomo with capability to calulculate pdfs
%cva_view_dots.m      ---  cva_slomo with capability to calulculate pdfs
%cvaf.m		      ---  like cva but with Fourier decomp (see cva_analysisf)
%demonstrate.m	      ---  demo 
%demonstrate3.m	      ---  demo
%fseries.m	      ---  n point, o order Fourier series
%pepsi.m	      ---  compare radial basis and polynomial techniques
%setpath.m            ---  setpath for the demos.





