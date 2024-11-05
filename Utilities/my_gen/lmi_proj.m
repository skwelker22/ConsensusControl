function [K,Ub]=lmi_proj(x0,A,B,RA,CE,toler,flag);
% function x=lmi_proj(x0,A,B,RA,CE,toler);
%  Projection via Convex optimization on the intersection
%  of ellipsoids (RA,CE) and half-spaces (A,B)
%

if nargin < 6; toler=0;flag=-1;end
if nargin < 7; flag=-1;end
if toler==0;toler=1.e-4;end

n=length(x0);[nx,me]=size(CE);mh=length(B);
Idn=eye(n,n);
K=x0;K0=hepprjn(x0,A,B,RA,CE,toler,1,flag);
err=sqrt(abs((K-K0)'*(K-K0)));
Aell=Idn*err*err*2;
U=1.e50;erro=1;kount=0;phi=0;kkount=0;indp=0;
[d2,dinf,D]=hseldist(K,A,B,RA,CE,1,flag);

while err>toler
  if d2>toler
      % 'constraint-iteration'
      indexx=find(D==dinf);index=indexx(1);
      iterty=-1;deepcut=0;
         if index<=mh
           h=A(index,:)';
         else
           h=2*RA(:,(index-mh-1)*n+1:(index-mh)*n)*(K-CE(:,index-mh));
         end
  else
      % 'objective-iteration'
      iterty=1; 
      h=2*(K-x0);phi=h'*h/4;h=h/sqrt(abs(h'*h));deepcut=0;
      erro=sqrt(abs(h'*Aell*h));
      if U>phi;U=phi;end
    end
[K,Aell]=lmi_upd(K,h,Aell,deepcut,flag);
kount=kount+1;kkount=kkount+1;
  [d2,dinf,D]=hseldist(K,A,B,RA,CE,1,flag);
  err=(d2)+erro;
end
Ub=sqrt(U);


