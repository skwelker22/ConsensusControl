function [L,R,Ub,x,k,UB,LB]=dia_con(A,i_rad,toler,x0,kkm);
% function [L,R,Ub,x]=dia_con(A,i_rad,toler,x0,kkm);
%   optimal conditioning of the matrix A using diagonal matrices
%   min cond(LAR); final R is scaled so that norm(LAR)=norm(A);


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
if length(x0)==n+m;x=x0;else;x=ones(n+m,1);end
COV=eye(n+m,n+m)*(i_rad^2);
Q=diag(x(1:m));P=diag(x(m+1:m+n));
error=toler+1;
k=0;kk=0;UBLB=1;
while error>toler
   Ab=A'*P*A;C=Ab-Q;
   [V,D]=eig(Q);dmin=min(real(diag(D)));index=max(find(real(diag(D))==dmin));
   v=V(:,index);v=v/norm(v);
   if dmin < 0
      g=zeros(m+n,1);iterty=-1;
      for i=1:m; g(i)=-v'*unit_vec(i,m)*unit_vec(i,m)'*v; end
      deepcut=-v'*Q*v;
   else
      [V,D]=eig(P);dpmin=min(real(diag(D)));index=max(find(real(diag(D))==dpmin));
      v=V(:,index);v=v/norm(v);
      if dpmin < 0
         g=zeros(m+n,1);iterty=-1;
         for i=1:n; g(i+m)=-v'*unit_vec(i,n)*unit_vec(i,n)'*v; end
         deepcut=-v'*P*v;
      else
         [V,D]=eig(C);dcmin=min(real(diag(D)));index=max(find(real(diag(D))==dcmin));
         v=V(:,index);v=v/norm(v);
         if dcmin < 0
            g=zeros(m+n,1);iterty=-1;
            for i=1:m; g(i)=v'*unit_vec(i,m)*unit_vec(i,m)'*v; end
            for i=1:n; g(i+m)=-v'*A'*unit_vec(i,n)*unit_vec(i,n)'*A*v; end
            deepcut=-v'*C*v;
         else
            [V,D]=eig(Ab,Q);dmax=max(real(diag(D)));index=max(find(real(diag(D))==dmax));
            v=V(:,index);v=v/norm(v);
            g=zeros(m+n,1);iterty=1;
            for i=1:m; g(i)=-dmax*v'*unit_vec(i,m)*unit_vec(i,m)'*v; end
            for i=1:n; g(i+m)=v'*A'*unit_vec(i,n)*unit_vec(i,n)'*A*v; end
            deepcut=0;
         end
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
      Q=diag(x(1:m));P=diag(x(m+1:m+n));
      error=max(error,UBLB);
      if any(x<0);error=max(error,toler+1);end
      if err_fl<0;error=0;end
      if kk>kkm
        disp([k,error,LB,UB]);kk=0;
      end
end
Ub=min(sqrt(UB));
R=sqrt(inv(Q));L=sqrt(P);
rat=norm(A)/norm(L*A*R);
R=R*rat;
