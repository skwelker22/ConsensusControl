function [J,DJ,DDJ]=feas_obj(x,A,B,RA,CE,indc);
%  function [J,DJ,DDJ]=feas_obj(x,A,B,RA,CE,indc);
%  P0-oblique projection of the vector x0 on the intersection
%  of ellipsoids RA=[R1,R2,...], CE=[c1,c2,...]
%  and half-spaces A,B

n=length(x);mh=length(B);[nx,me]=size(CE);mmh=length(A);mme=length(RA);
if mh==1 & mmh==1
  if A(1)==0 & B(1)==0, mh=0;end
end
if me==1 & mme==1
  if RA(1)==0 & CE(1)==0, me=0;end
end

  D=zeros(mh+me,1);k=0;
  DD=zeros(mh+me,n);
  DJ=0*x;
  DDJ=zeros(n,n);

  if mh>0
    zz=A*x-B;
    act_ind=find(zz>0);k=length(act_ind);
    D(1:k)=zz(act_ind);
       if indc(1)==1;DD(1:k,:)=A(act_ind,:);end
  end
  if me>0
    for i=1:me
      xperp=RA(:,(i-1)*n+1:i*n)*(x-CE(:,i));
      zzx=(abs((x-CE(:,i))'*xperp))-1;
        if zzx(1)>0;k=k+1;D(k)=zzx;
           if indc(1)==1;DD(k,:)=2*xperp';end
        end
    end
  end

if k>0
   J=D(1:k)'*D(1:k);
   if indc(1)==1
     DJ=2*(DD(1:k,:)'*D(1:k));
     DDJ=2*DD(1:k,:)'*DD(1:k,:);
   end
else
   J=0;
end