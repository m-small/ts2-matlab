function cva_movie(yp,ep,yt,pl,len)

% cva_movie(yp,ep,yt,pl,len)  draw some pretty pictures
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
if na<5
  len=pl;
  pl=[];
end;

figure(gcf);

miny=min([min(yt(:)),min(pl)]);
maxy=max([max(yt(:)),max(pl)]);
ps=length(yt(:,1));

global pl_pred_vect
ppv=pl_pred_vect;
ppv=-ppv;%(ppv<=0);


axis([0 ps+1 miny maxy]);
set(gca,'DrawMode','normal');

if na==5,
  if isempty(ep),
    i=1;
    pt=plot(yt(:,i),'k-.');hold on;pp=plot(yp(:,i));
    pv=plot(ppv+1,pl(i,:),'rp');  
    set(pt,'EraseMode','background');
    set(pp,'EraseMode','background');
    set(pv,'EraseMode','background');
    hold off;title(['Datum:',int2str(i)]);
    xlabel('data (black dot-dashed); cva-prediction (blue solid); radial basis prediction (red stars)');
    drawnow;
    axis([0 ps+1 miny maxy]);
  for i=1:len,
    set(pt,'YData',yt(:,i));
    set(pp,'YData',yp(:,i));
    set(pv,'YData',pl(i,:));
    set(get(gca,'Title'),'String',['Datum:',int2str(i)]);
    drawnow;
  end;
 else,
    i=1;
    pt=plot(yt(:,i),'k-.');hold on;pp=plot(yp(:,i));
    pph=plot(yp(:,i)+ep(:,i),'g-');
    ppl=plot(yp(:,i)-ep(:,i),'g-');
    pv=plot(ppv+1,pl(i,:),'rp');  
    set(pt,'EraseMode','background');
    set(pph,'EraseMode','background');
    set(ppl,'EraseMode','background');
    set(pp,'EraseMode','background');
    set(pv,'EraseMode','background');
    hold off;title(['Datum:',int2str(i)]);
    xlabel('data (black dot-dashed); cva-prediction (blue solid) and Std. Error (green solid); radial basis prediction (red stars)');
    drawnow;
    axis([0 ps+1 miny maxy]);
  for i=1:len,
    set(pt,'YData',yt(:,i));
    set(pp,'YData',yp(:,i));
    set(pph,'YData',yp(:,i)+ep(:,i));
    set(ppl,'YData',yp(:,i)-ep(:,i));
    set(pv,'YData',pl(i,:));
    set(get(gca,'Title'),'String',['Datum:',int2str(i)]);
    drawnow;
  end;
  end;
elseif na==4,
  i=1;
    pt=plot(yt(:,i),'k-.');hold on;pp=plot(yp(:,i));
    set(pt,'EraseMode','background');
    set(pp,'EraseMode','background');
    hold off;title(['Datum:',int2str(i)]);drawnow;
  for i=1:len,
    set(pt,'YData',yt(:,i));
    set(pp,'YData',yp(:,i));
    set(get(gca,'Title'),'String',['Datum:',int2str(i)]);
    drawnow;
  end;
end;

% End of cva_movie.m
