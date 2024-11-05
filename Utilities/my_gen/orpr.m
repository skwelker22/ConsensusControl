function [x,nitert]=orpr(x0,P0,A,B,RA,CE,fast,tol,flag,lam00);
%  function x=orpr(x0,P0,A,B,RA,CE,fast,tol,flag,lam00);
%  orthogonal projection of the vector x0 on the intersection
%  of ellipsoids RA=[R1,R2,...], CE=[c1,c2,...]
%  and half-spaces A,B
%  P0 defines the weighted distance x'P0x
%  optional arguments: tol=stopping tolerance (optim. tol=tol/100)
%  fast: >0 for hessian positivity check
%        =0 for adding 1e-5*In
%        <0 for no correction
%  lam00=initial and final point for constraint weight

flag=0;
if nargin<7;fast=1;end;if isempty(fast);fast=1;end
if nargin<8;tol=0;end;if isempty(tol);tol=0;end
if nargin<9;flag=0;end
if nargin<10;lam00=0;end
if tol(1)==0;tol=1.e-4;end
MITER=500;dp=tol*10;tloops=0;dtol=tol/100;

if length(lam00)==1 && lam00(1)==0;
   lam0=[1.e-2,1.e-8];
elseif length(lam00)==2
   lam0=lam00;
else
   lam0=[lam00,lam00*1e-6];
end

mh=length(B);mmh=length(A);
if mh==1 && mmh==1
   if A(1)==0 && B(1)==0, mh=0;end
end

%if mh >0
%   for i=1:mh
%     anor=sqrt(A(i,:)*A(i,:)');
%       if anor(1) > thresh
%         A(i,:)=A(i,:)/anor;B(i)=B(i)/anor;
%       else
%         A(i,:)=0*A(i,:);B(i)=0;
%          if flag(1) >1;disp('ignoring possibly singular constraint');end
%       end
%   end
%end

[d2,dinf,Dtem]=hseldist(x0,A,B,RA,CE,1,flag);
if dinf(1)<tol(1)
   x=x0;nitert=[0 0 0 0];
else
   xi=fgnewt(x0,A,(B-dp),RA*((1+dp)^2),CE,MITER,tol,fast,flag);    
   [d2,dinf,Dtem]=hseldist(xi,A,B,RA,CE,1,flag);
   lam1=lam0(1);x=xi;
	if flag>=1;disp('feasible point found');end
   if dinf(1)>0
      if flag(1) >0;disp('possibly empty constraint set');end
      nitert=[-1 -1 -1 -1];
      return
   end
   [x,niter1]=ognewt(xi,x0,P0,A,B,RA,CE,MITER,lam1,dtol,fast,flag);
	if flag>=1;disp('1st point found');end
   xi=x;xi1=x;
   lam2=lam1*.7;
   [x,niter2]=ognewt(xi,x0,P0,A,B,RA,CE,MITER,lam2,dtol,fast,flag);xi2=x;
	if flag>=1;disp('2nd point found');end
   niter=niter1+niter2;
   lam=lam2;tollam=tol*(lam1/lam2-1);dxn=10+tollam;indrun=0;

   while lam > lam0(2) && dxn(1) > tollam(1)
      indrun=indrun+1;tloops=tloops+1;
      lam=lam*(.02+.02^((indrun)/3));
      atemp=(xi2-xi1)/(lam2-lam1);dellam=lam-lam2;
      xi=x+atemp*dellam;
      [d2x,dinfx,Dtem]=hseldist(xi,A,B,RA,CE,1,flag);
      if dinfx(1)>0,
         while dinfx(1)>0,
            dellam=dellam*.8;lam=lam2+dellam;
            xi=xi2+atemp*dellam;
            [d2x,dinfx,Dtem]=hseldist(xi,A,B,RA,CE,1,flag);
         end
      elseif dinfx(1)<0
         xi=fgnewt(x0,A,(B-dp),RA*((1+dp)^2),CE,MITER,tol,fast,flag);    
         [d2x,dinfx,Dtem]=hseldist(xi,A,B,RA,CE,1,flag);
      end
      [x,niter3]=ognewt(xi,x0,P0,A,B,RA,CE,MITER,lam,dtol,fast,flag);
      niter=niter+niter3;
      xi1=xi2;xi2=x;lam1=lam2;lam2=lam;
      dxn=sqrt(abs((xi2-xi1)'*(xi2-xi1)));tollam=tol*(lam1-lam2)/lam2;
	if flag>=1
		disp('dx/tol_lam      n_iter     lam')
		disp([dxn/tollam, niter lam])
	end
      if abs(lam-lam1)<tol*lam1;
         dxn=0;
         disp(['orpr: forced exit; lam = ',num2str(lam)])
      end
   end
   nitert=[niter1,niter2,tloops,niter];
end

