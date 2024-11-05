function [lspro,index]=lsproj(th,thmax,thmin,P,PI0,projtol);
% usage: [lspro,index]=lsproj(th,thmax,thmin,P,PI,projtol);
% projection on half spaces for least-squares adaptation
% index=1 if projection is active

MAXIT=500;
if nargin<4
   P=0;PI0=0;projtol=0;
elseif nargin<5
   PI0=0;projtol=0;
elseif nargin < 6
   projtol=0;
end
if projtol(1)==0;  projtol=1.e-6;end
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

index=0;lspro=th;A=[It;-It];B=[tmax;-tmin];
er1=A*th-B;err=max(max(er1),0);
kiter=0;

while err>projtol(1,1)
   lspr=lspro;kiter=kiter+1;
   for ipar=1:2*n
      aa=A((ipar),:);vp=aa*PI;
      num=min(B((ipar)),aa*lspro)-aa*lspro;
      lspro=lspro+num/(vp*aa')*vp';
   end
   er1=A*lspro-B;err=max(max(er1),0);
   if err(1,1)>projtol(1,1)
      index=1;
   end
   if kiter>MAXIT;
      disp(['lsproj: max iterations exceeded.  err = ',num2str(err)]);
      th
      PI
      err=0;
   end
end
