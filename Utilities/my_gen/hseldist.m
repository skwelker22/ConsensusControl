function [d2,dinf,D]=hseldist(x,A,B,RA,CE,dista,flag);
%
% USAGE: [d2,dinf,D]=hseldist(x,A,B,RA,CE,dista);
%
% Computes the distance of x from half-spaces or ellipsoids
% d2=2-norm of distances, dinf=max distance 
% dista=0 for distance, 1 for equation error
% flag<=1 for silent run

n=length(x);mh=length(B);[nx,me]=size(CE);mmh=length(A);mme=length(RA);
if mh==1 & mmh==1
  if A(1)==0 & B(1)==0, mh=0;end
end
if me==1 & mme==1
  if RA(1)==0 & CE(1)==0, me=0;end
end
D=zeros(mh+me,1);d2=0;dinf=0;

if mh>0
   z=A*x-B;
   act_ind=find(z>0);k=length(act_ind);
   if dista(1) ~=0
      D(act_ind)=z(act_ind);

   else
      for i=1:k
        ii=act_ind(i);
        ztem2=sqrt(abs(A(ii,:)*A(ii,:)'));
         if ztem2(1)<1.e-12
           if flag(1) >1;disp('ignoring possibly singular constraint');end
         else
           D(ii)=z(ii)/ztem2;
         end
      end
   end
end
if me>0
   for i=1:me
     xperp=RA(:,(i-1)*n+1:i*n)*(x-CE(:,i));
     xce=sqrt(abs((x-CE(:,i))'*xperp))-1;
     if D(mh+i)< xce(1);D(mh+i)=xce;end
   end
end
d2=sqrt(abs(D'*D));
dinf=max(D);

%for i=1:mh+me;
%  if dinf<D(i);dinf=D(i);end
%end

