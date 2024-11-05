function [L,Ub,x,k,UB,LB]=dia_conl(A,i_rad,toler,x0,kkm);
% function [L,Ub,x]=dia_conl(A,i_rad,toler,x0,kkm);
%   optimal conditioning of the matrix A using diagonal matrices
%   min cond(LA); final L is scaled so that norm(LA)=norm(A);

if nargin<2;i_rad=0;end
if nargin<3;toler=0;end
if nargin<4;x0=0;end
if nargin<5;kkm=49;end
if isempty(i_rad) | i_rad==0;i_rad=100;end
if isempty(toler) | toler==0;toler=0.1;end


L=[];R=[];dmin=0;
[n,m]=size(A);
if n <m;disp('dia_con: matrix A must have more rows than columns');return;end
dmax=cond(A);LB=[];UB=[];
if length(x0)==n;x=x0;else;x=ones(n,1);end
COV=eye(n,n)*(i_rad^2);Q=eye(m,m);P=diag(x);
error=toler+1;
k=0;kk=0;UBLB=1;
while error>toler
   Ab=A'*P*A;C=Ab-Q;
   [V,D]=eig(P);dmin=min(real(diag(D)));index=max(find(real(diag(D))==dmin));
   v=V(:,index);v=v/norm(v);
   if dmin < 0
      g=zeros(n,1);iterty=-1;
      for i=1:n; g(i)=-v'*unit_vec(i,n)*unit_vec(i,n)'*v; end
      deepcut=-v'*P*v;
   else
      [V,D]=eig(C);dcmin=min(real(diag(D)));index=max(find(real(diag(D))==dcmin));
      v=V(:,index);v=v/norm(v);
      if dcmin < 0
         g=zeros(n,1);iterty=-1;
         for i=1:n; g(i)=-v'*A'*unit_vec(i,n)*unit_vec(i,n)'*A*v; end
         deepcut=-v'*C*v;
      else
         [V,D]=eig(Ab);dmax=max(real(diag(D)));index=max(find(real(diag(D))==dmax));
         v=V(:,index);v=v/norm(v);
         g=zeros(n,1);iterty=1;
         for i=1:n; g(i)=v'*A'*unit_vec(i,n)*unit_vec(i,n)'*A*v; end
         deepcut=0;
      end
   end

%      g=g/norm(g);
      k=k+1;kk=kk+1;
      if iterty>0
         error=sqrt(abs(g'*COV*g));
         if length(UB)>=1
            LB=max(LB,real(dmax)-error);UB=min(UB,real(dmax));
         else
            LB=real(dmax)-error;UB=real(dmax);
         end
         UBLB=UB-LB;
%          LB=[LB;dmax-error];UB=[UB;dmax];
      end
      [x,COV,err_fl]=LMI_upd(x,g,COV,deepcut);x=real(x);COV=real((COV+COV')/2);
      P=diag(x(1:n));
      error=max(error,UBLB);
      if any(x<0);error=max(error,toler+1);end
      if err_fl<0;error=0;end
      if kk>kkm
        disp([k,error,LB,UB]);kk=0;
      end
      if kk>kkm
        disp([k,error,LB,UB]);kk=0;
      end
end
Ub=min(sqrt(UB));
L=sqrt((P));
rat=norm(A)/norm(L*A);
L=L*rat;
