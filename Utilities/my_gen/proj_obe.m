function [J,DJ,DDJ]=proj_obe(x,x0,P0,A,B,RA,CE,lam,indc);
%  function [J,DJ,DDJ]=proj_obe(x,x0,P0,A,B,RA,CE,lam,indc);
%  P0-oblique projection of the vector x0 on the intersection
%  of ellipsoids RA=[R1,R2,...], CE=[c1,c2,...]
%  and half-spaces A,B

if lam(1)<=0;lam=1.e-15;end
if length(P0)==1 & P0(1,1)==0;P0=1;end

n=length(x);mh=length(B);[nx,me]=size(CE);mmh=length(A);mme=length(RA);
if mh==1 & mmh==1
   if A(1)==0 & B(1)==0, mh=0;end
end
if me==1 & mme==1
   if RA(1)==0 & CE(1)==0, me=0;end
end

ii=[1:n];iii=ii+n*(ii-1);
Idn=zeros(n,n);Idn(iii)=Idn(iii)+1;
D=zeros(mh+me,1);d2=0;DD=zeros(mh+me,n);ONV=D+1;DDX=DD;
DJ=0*x;DDJ=Idn;DL=D;

if mh>0
   D(1:mh)=A*x-B;
   if indc(1)==1;DD(1:mh,:)=A;end
end
if me>0
   for i=1:me
      xperp=RA(:,(i-1)*n+1:i*n)*(x-CE(:,i));
      D(mh+i)=(abs((x-CE(:,i))'*xperp))-1;
      if indc(1)==1;DD(mh+i,:)=2*xperp';end
   end
end  

DL=0*D;
ind_D=find(D>=0);
if isempty(ind_D)
   J=1.e51;
else
   DL(ind_D)=D(ind_D).*D(ind_D);
   J=(x-x0)'*P0*(x-x0)+lam*ONV'*DL;
   if indc(1)==1
      DJ=2*P0*(x-x0)+2*lam*(DD(ind_D,:)'*(D(ind_D,:)));
      DDJ=2*P0*Idn+2*lam*DD(ind_D,:)'*DD(ind_D,:);
      if me>0
         for i=1:me
            DDJ=DDJ+2*lam*RA(:,(i-1)*n+1:i*n);
         end
      end   
   end   
end
