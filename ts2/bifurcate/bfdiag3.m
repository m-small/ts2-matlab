function [yt,y,ps,ts,pts]=bfdiag3(y,mu,h,n,res);

%function [yt,yp,p,t,pt]=bfdiag3(y,mu,h,n,res)
%
% same as bfdiag3, but works with really looong predicted runs.
%
% do calculations and plot a bifurcation diagram
%
% make predictions h points into the future, and compute n peak/troughs for each future prediction
% make at most k future predictions. res is the upper bound on how
% many i.c. to use
%
% yt is the test data time series, yp are the k future predictions, p and t are the n peak and 
% trough values for each of the k future predictions, and pt are the datum number (mapped to yt)
% for which the long term behaviour is computed.
%
% If mu is empty the use constant increase over the range
% provided. If y is empty then use rb_x for both y (embedded) and
% mu.
%
% defualts:
% h=1000
% n=10
% k=500
%
% M. Small 
% Created: 17/10/99
% Updated: 19/10/99

na=nargin;
if na<5,
  res=500;
end;
if na<4,
   n=10;
   if na<3
      h=1000;
      if na<2,
         mu=[];
         if na<1,
           y=[];
	 end;
      end;     
   end;
end;

%block size, bb, number of iterations to do simulataneously
bb=100;

rb_get_globals;
v=rb_embed;
v=v(~isnan(v(:)));
v=0:max(v);


[dy,ly]=size(y);
rescale=1;
if ly==1,
   y=y';
   [dy,ly]=size(y);
elseif (dy==0 & ly==0),
  y=rb_x;
  [dy,ly]=size(y);
   mu=y(dy,:);
   rescale=0;
end;

if isempty(mu),
   mu=1:ly;
end;

if rescale,
  mu=bfmu(mu);
end;

%sort out the initial conditions
if dy==1,
   yt=y;
   y=embed(y,v);
elseif dy~=max(v),
   yt=y(1,:);
   disp('WARNING : Dimensions don''t agree');
end;

[dy,ly]=size(y);
if length(mu)<ly,
   disp('WARNING : size of mu and y don''t agree');
end;

%shorten y and mu, so that we don't do too much work
if ly>res,
   ss=floor(ly./res);
   y=y(:,1:ss:end);
   mu=mu(1:ss:end);
   [dy,ly]=size(y);
end;


%do the freeruns & plots 
subplot(311);
plot(yt);
title('Data');axis tight
subplot(312);hold on;title('Peaks (green) and Troughs (red)')
subplot(313);hold on;title('Amplitudes');
ps=nan*ones(ly,n+1);ts=ps;pts=[];
y=y(1:dy,:);
for i=1:bb:ly
  ir=i+(1:bb)-1;
  ir=ir(ir<=ly);
  yy=y(:,ir);
  yy=bfiterate(yy,mu(ir),h);
  for j=1:(ir(end)-i+1),
   [p,pt,t,tt]=peaktrough(yy(:,j),10,0);
%[p,pt,t,tt]=extremum(y(:,i),[],[],[],1,1);
   lp=length(p);
   if lp>0,
      lp=max(lp-n,1);
      p=p(lp:end);
      t=t(lp:end);
      lp=length(p);
      ps(i+j-1,1:lp)=p;
      ts(i+j-1,1:lp)=t;
      subplot(312);
      plot((i+j-1)*ss*ones(size(p)),p,'g.',(i+j-1)*ss*ones(size(p)),t,'r.');
      subplot(313);plot((i+j-1)*ss*ones(size(p)),p-t,'c.');
   end;
   xlabel(['done : ',int2str((i+j-1)),'/',int2str(ly)]);
   drawnow;
   end;
end;
xlabel('');
subplot(312);axis tight;box on;
subplot(313);axis tight;box on;
pts=(1:ly)*ss;
