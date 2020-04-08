function [z,zx,zy]=pdf2d(x,bw,res,xlim,ylim)

%create 2d pdf from the discrete observations in the 2-by-n vector x
%bw is the optional bandwidth parameter for kernel smoothing of the result
%res is the resolution of the output
%

if nargin<5,
    ylim=[];
end;
if nargin<4,
    xlim=[];
end;

if nargin<3,
    res=100;
end;
if nargin<2,
    bw=1;
end;

[dx,nx]=size(x);
if dx~=2 & nx==2,
    x=x';
    [dx,nx]=size(x);
end;

if isempty(xlim),
    minx=min(x(1,:));
    maxx=max(x(1,:));
else
    minx=min(xlim);
    maxx=max(xlim);
end;
if isempty(ylim),
miny=min(x(2,:));
maxy=max(x(2,:));
else
    miny=min(ylim);
    maxy=max(ylim);
end;
stpx=(maxx-minx)/res;
stpy=(maxy-miny)/res;

[zx,zy]=meshgrid(minx:stpx:maxx,miny:stpy:maxy);
z=0*zx;

for i=1:nx,
    z=z+exp( -(  ((zx-x(1,i))/(stpx*bw)).^2 + ((zy-x(2,i))/(stpy*bw)).^2 )/2 );
end;
    
z=z./nx;