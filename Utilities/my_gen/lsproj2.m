function [lsproj,index]=lsproj2(th,thmax,thmin,P,PI0,projtol);
% usage: [lsproj,index]=lsproj2(th,thmax,thmin,P,PI,projtol);
% projection on half spaces for least-squares adaptation
% index=1 if projection is active

nma=length(thmax);nmi=length(thmin);n=length(th);
It=eye(n,n);

if length(P)==1 & P(1,1)==0;  P=1;end
if length(PI0)==1 & PI0(1,1)==0;  PI=P\It;else PI=PI0;end
if nma==n
   tmax=thmax;
elseif nma==1
   tmax=zeros(n,1)+thmax;
else
   tmax=zeros(n,1)+1.e50;
end
if nmi==n
   tmin=thmin;
elseif nmi==1
   tmin=zeros(n,1)+thmin;
else
   tmin=zeros(n,1)-1.e50;
end

index=0;lsproj=th;A=[It;-It];B=[tmax;-tmin];
er1=A*th-B;err=max(max(er1),0);

if err>projtol(1,1)
   index=1;
   %   PX=P/norm(P);
%   p1=sqrtm(PI);p2=sqrtm(P);
%   A=A*p1;th=p2*th;
%   lsproj=orpr(th,1,A,B,[],[],1,projtol,0);
%   lsproj=p1*lsproj;
lsproj=min(th,tmax);
lsproj=max(lsproj,tmin);
end
