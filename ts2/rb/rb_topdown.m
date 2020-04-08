function [base,best_a,best_deltas,best_basis,e,mdl]= rb_topdown(X,y,v,func,base,penalty,nc,maxcount,wobble,meth)

% function [base,lambdas,deltas,basis,err,mdl]= rb_topdown(X,y,v,func,base,penalty,nc,maxcount,wobble,meth)
%
% calculate possible L2-norm best radial basis model of y from cols of X
% using sensitivity analysis and MDL
% base is initial candidate basis functions
% maxcount= max. # of additions/substitutions of basis funcs that do not decrease MDL
% wobble, if wobble then perform optimisation (Nelder-Mead Simplex descent) on each new basis candidate.
% penalty is a string specifying the penalty equation to eval e.g.
%   Akaike = 'dim_x*log(mss)+2*k'                                               --- Akaike criterion
%   Schwarz = 'dim_x*log(mss)+k*log(dim_x)',                                    --- Schwarz criterion
%   SRissanen = 'dl1(mss,unscaled_a,unscaled_deltas,dim_x)';                     --- pseudo-linear DL
%   SRissanen2= 'dl2(mss,unscaled_a,unscaled_deltas,dim_x,base,basic,offset,v)'; --- approx DL
%   Rissanen = 'dl1(mss,a,deltas,dim_x)';                     --- pseudo-linear DL (unscaled params)
%   Rissanen2= 'dl2(mss,a,deltas,dim_x,base,basic,offset,v)'; --- approx DL (unscaled params)
%
%
% base is a list of basis functions selected
% err are residual errors
% mdl the normalized description length of model
%
% best_a is a nb-by-2 matrix, the first column is the (scaled) best weights and the second column is
% the scaling factor. Similarly best_deltas are scaled by the same factor.
% This prevents (i.e. reduces) numerical sensitivity.
%
% Best model is y = rb_eval(X)+err = rb_Phi(X,base,v)*base.lambda+err;
%
% Defaults:
%   penalty = Rissanen
%   maxcount= 5
%   basic=[];
%   wobble=0;
%   meth='clr'
%
% trace= 0 no display & no info
% trace= 1 display & basic info
% trace< 0 info as above but no display
%
% M. Small 
% Created: 7/3/99
% Updated: 8/4/05

% penalty functions
Akaike = 'dim_x*log(mss)+2*k';                %AIC
Schwarz = 'dim_x*log(mss)+k*log(dim_x)';      %SIC/BIC
Rissanen = 'dl1(mss,a,deltas,dim_x)';  %pseudo-linear DL
%i.e. penalise for description of wieghts only
Rissanen2= 'dl2(mss,a,deltas,dim_x,base,basic,offset,v)'; %approx DL
%Why not unscaled_a and unscaled_deltas!?
%assume penalty for all significant parameters are equal and equal to the penalty for
%the weights (this assumption is not necessarily well founded --- just convenient)

%Exact_Rissanen = 'description_length2(mss,unscaled_a,unscaled_deltas,dim_x)';
%Really_Exact_Rissanen = 'description_length2(X,y,best_basis,unscaled_a,unscaled_deltas)';

a=[];r=[];

y= y(:);
[nx,dim_x]= size(X);
[dim_y,ny]= size(y);

if dim_y~=dim_x
   error('dim_y~=dim_x');
end
if ny~=1
   error('ny~=1');
end
if nargin<3,
   v=[-1 0 1 2];
end;
if nargin<4,
   func{1}='gaussian';
end;
if nargin<5
   base=rb_nullbase;
end
if nargin<6
   penalty= Rissanen;
end
if nargin<7,
   nc=300;
end;
if nargin<8
   maxcount= 5;
end
if nargin<9,
   wobble=0;
end;
if nargin<10,
  meth='clr';
end;

offset=0;
if any(meth=='c'),
  offset=offset+1;
end;
if any(meth=='l'),
  offset=offset+nx;
end;

%initial model
if isempty(base),
   base=rb_nullbase;
end;
obase= base; 
no=length(base.func);
newold=[];
deleted=[];
base=rb_nullbase;
phi=rb_Phi(X,obase,v,func,meth);
[phi,scale]= normalize(phi);                % necessary to avoid ill-conditioning - definitely
basic=[];
unscaled_deltas=[];
expansion=0;

%should we bother with precisions?
if strcmp(penalty,Schwarz) | strcmp(penalty,Akaike),
   needdeltas=0;
else,
   needdeltas=1;
end;

%compute screensize and figure position
figpos=get(0,'ScreenSize');
figpos=floor([20 3*figpos(4)/5-20 min(figpos(3),1800)-40 figpos(4)/3]);

% display
if trace>0
   fig1 = figure;
   disp(sprintf('Displaying on figure %d',fig1.Number));
   set(fig1,'Position',figpos);
   subplot(211);
   handle1=plot(y,'b');hold on;
   plot(y,'k');hold off;
   set(gca,'XTickLabel','');
   drawnow;
   subplot(212);
   handle2=plot(y,'b'); 			% the current error
end

%information
if trace>=1
   disp('   No.     RMS        DL    :  basis...');
end;


% normalize X and retain scale factor

% set up X=QR and initial deltas
k= length(basic);
nl= length(base.strategy);
if k==0
   q= [];
   deltas= [];
else
   [q,r]= qr(phi(:,basic-offset),0);
   if ~isempty(base.deltas)
      deltas= base.deltas .* scale(basic)'; 	% scale supplied deltas
   end
end



mdl= Inf; 												% minimum description length so far
next= []; 												% next vector to add to basis
count= 0; 												% number of failed expansions

pmss=inf;                                    % previous mss: trap to prevent switch of 
%                                              basis functions leading to increase in
%                                              mss (due to nonorthogonality of candidates)


while k<ny
   
   %k= length(basic); 			               % current model size
   improvement= 0; 									% no improvements yet
   all_basis=basic;
   
   j= k; 			   	         					% if j> k, then it is time
   while j <= k 										% to expand model
      
      % evaluate model
      phi=rb_Phi(X,base,v,func,meth);
      [phi,scale]=normalize(phi);
      
      if isempty(basic)
         e= y; 											% current error
         mss= e'*e/dim_y; 								% mean square error
         S= [];
      else
         a= r\(q'*y); 									% parameters
         e= y - phi(:,basic)*a; 						% current error
         mss= e'*e/dim_y; 								% mean square error
         S= r'*r/mss; 									% second derivative of log-likelihood
      end
      
      %check for swap leading to increase in mss,
      if pmss<mss
      %   disp('mss rise - caught');
         %belated expansion of basis
         expansion=1;
         nl=nl+1;
         base=pbase;basic=pbasic;
         q=pq;r=pr;
         mss=pmss;
         newold=pnewold;
         deleted=[];
         %expand
         k=k+1;
         if deltas,
            deltas=[deltas;mean(deltas)];
         end;
         %recalculate phi
         phi=rb_Phi(X,base,v,func,meth);
         [phi,scale]=normalize(phi);
         %recalculate mss, q, r and S
         a= r\(q'*y); 									
         e= y - phi(:,basic)*a;
         mss= e'*e/dim_y;
         S= r'*r/mss;
      end;
      %reset variables for next time
      pbase=[];pbasic=[];
      pq=[];pr=[];
      pnewold=[];
      pmss=inf;
     
      
      if ~isempty(S) & needdeltas
         deltas= find_deltas(S,deltas); 			% precisions of 
      end
      
      % current model's description length
      unscaled_a= a ./ scale(basic)';
      if needdeltas,
         unscaled_deltas= deltas ./ scale(basic)';   
      end;
      descr_length = eval(penalty);
      
      % Has model reduced description length?
      improvement=0;
      if mdl > descr_length
         mdl= descr_length; 							% keep best model info
         best_basis= basic;
         best_base=base;
         best_e= e;
         best_rms= sqrt(mss);
         best_a = [a scale(basic)'];% to avoid numerical error store normalized values
                                                % was: best_a=unscaled_a;
         best_deltas = deltas; 						% was: unscaled_deltas;
         improvement= 1; 								% extension has yielded improvement
      end
      
      
      
      
      % information
      if improvement,
         message0='*';
      else,
         message0=' ';
      end;
      message1= sprintf(' %3d %10.4G %10.4G :',k,sqrt(mss),descr_length);
      message2= [' ',rb_string(basic,base,offset,func,newold,deleted)];
      if trace>=1
         disp([message0 message1 message2]);
      end;
      
      % try improving basis by extend and delete  
      % extra candidates for extend
      if any(meth=='n'),
          %if its a neural net choose this way,
          nbase=rb_choose_net(X,[],func,v,nc);
      else,
          %otherwise do this
          nbase=rb_choose_new(X,[],func,v,nc);
      end;
      %now optimise radii?
%      nbase=rb_stretch_r(nbase,X,v,func,e,meth);
      %or everything?
%      nbase=rb_stretch_all(nbase,X,v,func,e,meth);
      
      %concatenate new and existing basis functions
      cbase=rb_catbase(base,obase,nbase);
      phi=rb_Phi(X,cbase,v,func,meth);
      [phi,scale]=normalize(phi);

      % first extend
      mu = -phi'*e;     								% calculate sensitivity mu of vectors
      [s,i] = max(abs(mu)); 							% find most sensitive
      basic = [ basic, i ]; 							% move it into basis
      if i>offset,             						% is new entrant nonlinear?
         nl=nl+1;											% if so, increase the count and 
         basic(k+1)=nl+offset;						% add it to the basis
         newold(k+1)=((i-offset)<=no);
         newbase=rb_assbase(cbase,i-offset);
         if wobble>0,       
            %attempt to refine the selected basis function
            newnewbase=newbase;
            [newbase,thisphi,thisscale,improved]=rb_wobble(newbase,X,v,func,e,meth,wobble);
           % if abs(phi(:,i)'*e)>abs(thisphi'*e),
	   %   disp('WARNING: wobbling failed!');
	   %   disp(['was=',num2str(-phi(:,i)'*e)]);
	   %   disp([' is=',num2str(-thisphi'*e)]);
	   % end;
	    scale(i)=thisscale;thisscale=[];
            newold(k+1)=newold(k+1)*improved;
         else,
            thisphi=phi(:,i);
         end;
         base=rb_catbase(base,newbase); %%i-offset->base(j)-offset
      else,
         thisphi=phi(:,i);
      end;
      [q1,r1] = qrappend(q,r,thisphi); 			% extend basis
      a1 = r1\(q1'*y); 								% new model parameters
      thisphi=[];
      
      % now try a delete
      previous=expansion;
      [s,j]= min(abs(a1)); 							% poorest vector
      if j==k+1,   										% delete vector just added?
         %     disp('last is worst');
         expansion=1;
         q= q1; 											% yes, then extend basis: so
         r= r1; 											% keep factorization and
         if deltas
            deltas= [deltas;mean(deltas)]; 		% add a likely precision
         end
         k=k+1;
         deleted=[];
      else,
         expansion=0;
         %    disp('last is NOT worst');
         if basic(j)>offset,
            %     disp('worst is nonlinear');
            nl=nl-1;
            pbase=base;pbasic=basic;
            pq=q1;pr=r1;pmss=mss;pnewold=newold;
            base=rb_delbase(base,basic(j)-offset); % delete from base set
            basic(basic>basic(j))=basic(basic>basic(j))-1; 		% adjust other indicies
            newold(j)=[];
         end;
         basic(j)=[];
         deleted=j;
         [q,r]= qrdelete(q1,r1,j); 					% and QR factorization
         q= q(:,1:k); 		       					% remove stuff that
         r= r(1:k,1:k); 			   				% qrdelete doesn't remo
      end
      
      % count is an ad hoc method to decide when to stop:
      % --   Finding a new mdl resets count
      % --   each time the model is extended without improvement of mdl
      % --   then count increases
      
      if any(basic>offset),
         % keyboard;
      end;
      
      % Has model reduced description length?
      if improvement
         count= 0;
      elseif previous %expansion
         count= count+1;
      end;
      %disp([previous,improvement,count,maxcount]);
      % extending model to no avail? Then quit
      if count >= maxcount
         break
      end
      
      % Global optimal solution?
      if  mss < eps,
         disp('Global optimal solution, or near enough');
         break;
      end;
      
      
      % nice display
      if trace>0 & improvement
         figure(fig1); subplot(212);
         set(handle1,'ydata',best_e);
         set(handle2,'ydata',best_e);
         message1= sprintf('%d basis functions, RMS= %G, MDL= %G', ...
            length(best_basis),best_rms,mdl);
         title(message1);
         drawnow;
      end
      
   end;
end;

if trace
   message1= sprintf(' %d basis functions, RMS= %G, MDL= %G', ...
      length(best_basis),best_rms,mdl);
   % message2= sprintf('%3d ',best_basis);
   disp(['Best compression is possibly at ' message1]);
   % disp(message2);
end


phi=rb_Phi(X,best_base,v,func,meth);
[phi,scale]=normalize(phi); %necessary to avoid ill-conditioning ???

if isempty(best_basis),
  e=y;
  p=e-y;
else
    p=phi(:,best_basis)*best_a(:,1);%.* scale(best_basis)');
    e= p-y; 				
end;

lambda = best_a;
precision = best_deltas;
base=best_base;

if trace>0
   clf, 
   plotcols(e','b:',y','k-',p','r--');
   title(message1);
   drawnow;
end
if rms(e')>rms(y'),
   keyboard;
end;

if 0 & ~needdeltas,
   best_deltas=e;
   e=mdl;
   mdl=[];
end;
