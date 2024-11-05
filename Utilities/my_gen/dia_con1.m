function [L,R,Ub,x,k,UB,LB]=dia_con1(A,x0,i_rad,toler);
% function [L,R,Ub,x]=dia_con1(A,x0,i_rad,toler);
%   optimal conditioning of the matrix A using diagonal matrices
%   min cond(LAR)

L=[];R=[];dmin=0;
[n,m]=size(A);
if n <m;disp('dia_con: matrix A must have more rows than columns');return;end
dmax=cond(A);LB=1;UB=dmax^2;
if length(x0)==n+m-1;x=x0;else;x=ones(n+m-1,1);end
COV=eye(n+m-1,n+m-1)*(i_rad^2);
Q=diag(x(1:m));P=diag([x(m+1:m+n-1);1]);
error=toler+1;
k=0;
while error>toler
   Ab=A'*P*A;C=Ab-Q;
   [V,D]=eig(Q);dmin=min(real(diag(D)));index=max(find(real(diag(D))==dmin));
   v=V(:,index);v=v/norm(v);
   if dmin < 0
      g=zeros(m+n-1,1);iterty=-1;
      for i=1:m; g(i)=-v'*unit_vec(i,m)*unit_vec(i,m)'*v; end
      deepcut=-real(v'*Q*v);
   else
      [V,D]=eig(P);dpmin=min(real(diag(D)));index=max(find(real(diag(D))==dpmin));
      v=V(:,index);v=v/norm(v);
      if dpmin < 0
         g=zeros(m+n-1,1);iterty=-1;
         for i=1:n-1; g(i+m)=-v'*unit_vec(i,n)*unit_vec(i,n)'*v; end
         deepcut=-real(v'*P*v);
      else
         [V,D]=eig(C);dcmin=min(real(diag(D)));index=max(find(real(diag(D))==dcmin));
         v=V(:,index);v=v/norm(v);
         if dcmin < 0
            g=zeros(m+n-1,1);iterty=-1;
            for i=1:m; g(i)=v'*unit_vec(i,m)*unit_vec(i,m)'*v; end
            for i=1:n-1; g(i+m)=-v'*A'*unit_vec(i,n)*unit_vec(i,n)'*A*v; end
            deepcut=-real(v'*C*v);
         else
            [V,D]=eig(Ab,Q);dmax=max(real(diag(D)));index=max(find(real(diag(D))==dmax));
            v=V(:,index);v=v/norm(v);
            g=zeros(m+n-1,1);iterty=1;
            for i=1:m; g(i)=-dmax*v'*unit_vec(i,m)*unit_vec(i,m)'*v; end
            for i=1:n-1; g(i+m)=v'*A'*unit_vec(i,n)*unit_vec(i,n)'*A*v; end
            deepcut=0;
         end
      end
   end

%      g=g/norm(g);
      k=k+1;
      if iterty>0
         error=sqrt(abs(g'*COV*g));
%         LB=max(LB,(dmax)-error);UB=min(UB,(dmax));
          LB=[LB;dmax-error];UB=[UB;dmax];
      end
      [x,COV]=LMI_upd(x,g,COV,deepcut);x=real(x);COV=real((COV+COV')/2);
      Q=diag(x(1:m));P=diag([x(m+1:m+n-1);1]);
%      error=max(error,UB-LB);
      if any(x<0);error=max(error,toler+1);end
     disp([k,error,iterty])
end
Ub=min(sqrt(UB));
R=chol(inv(Q));L=chol(P);
