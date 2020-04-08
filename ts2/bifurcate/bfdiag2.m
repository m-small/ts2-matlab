function [yt,y,ps,ts,pts]=bfdiag2(y,mu,h,n,res);

%function [yt,yp,p,t,pt]=bfdiag2(y,mu,h,n,res)
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
% Same as bfdiag, except that the i.c. are chosen randomly over the 
% embedding space (except for the bifurcation parameter
% coordinate). 
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


%do the freeruns
tn=min(h,200); %200 points at a time, just to save some memory.
while tn>0,
   y=y(1:dy,floor(rand(1,ly)*ly+1));
   y=bfiterate(y,mu,tn);
   h=h-tn;
   tn=min(h,200);
end;

%do the plots
subplot(311);
plot(yt);
title('Data');axis tight
subplot(312);hold on;title('Peaks (green) and Troughs (red)')
subplot(313);hold on;title('Amplitudes');
ps=nan*ones(ly,n+1);ts=ps;pts=[];
for i=1:ly
   [p,pt,t,tt]=peaktrough(y(:,i),10);
%[p,pt,t,tt]=extremum(y(:,i),[],[],[],1,1);
   lp=length(p);
   if lp>0,
      lp=max(lp-n,1);
      p=p(lp:end);
      t=t(lp:end);
      lp=length(p);
      ps(i,1:lp)=p;
      ts(i,1:lp)=t;
      subplot(312);plot(i*ss*ones(size(p)),p,'g.',i*ss*ones(size(p)),t,'r.');
      subplot(313);plot(i*ss*ones(size(p)),p-t,'c.');
   end;
   xlabel(['done : ',int2str(i),'/',int2str(ly)]);
   drawnow;
end;
xlabel('');
subplot(312);axis tight;box on;
subplot(313);axis tight;box on;
pts=(1:ly)*ss;
