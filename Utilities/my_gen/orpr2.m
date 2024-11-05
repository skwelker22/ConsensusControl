function [x,nitert]=orpr2(x0,P0,A,B,RA,CE,fast,tol,flag,lam00);
%  function x=orpr2(x0,P0,A,B,RA,CE,fast,tol,flag,lam00);
%  orthogonal projection of the vector x0 on the intersection
%  of ellipsoids RA=[R1,R2,...], CE=[c1,c2,...]
%  and half-spaces A,B
%  optional arguments: tol=stopping tolerance (optim. tol=tol/100)
%  fast: >0 for hessian positivity check
%        =0 for adding 1e-5*In
%        <0 for no correction
%  lam00=initial and final point for constraint weight

if nargin<7;fast=-1;end;if isempty(fast);fast=-1;end
if nargin<8;tol=0;end;if isempty(tol);tol=0;end
if nargin<9;flag=0;end
if nargin<10;lam00=0;end
if tol(1)==0;tol=1.e-4;end
MITER=200;tloops=0;dtol=tol/100;
BLR=0.6;lam_fact=.6;dinf_thr=1;

if length(lam00)==1 
   if lam00(1)==0;
      lam0=[1.e-1,1.e-12,lam_fact,dinf_thr];
   else
      lam0=[lam00,lam00*1e-11,lam_fact,dinf_thr];
   end
elseif length(lam00)==2
   lam0=[lam00,lam_fact,dinf_thr];
elseif length(lam00)==3
   lam0=[lam00,dinf_thr];
else
   lam0=lam00;
end
lam_fact=lam0(3);

mh=length(B);mmh=length(A);
if mh==1 && mmh==1
   if A(1)==0 && B(1)==0, mh=0;end
end

[d2,dinf,Dtem]=hseldist(x0,A,B,RA,CE,1,flag);
if dinf(1)<tol(1)
   x=x0;nitert=[0 0 0 0];
else
   lam1=lam0(1);
   BL=B+dinf*BLR*lam1;
   [x,niter1]=ognewte(x0,x0,P0,A,BL,RA,CE,MITER,1/lam1,dtol,fast,flag);
   if flag>=1;disp('1st point found');end
   xi=x;xi1=x;
   lam2=lam1*.8;
   [d2,dinf,Dtem]=hseldist(xi,A,B,RA,CE,1,flag);
   BL=B+dinf*BLR*lam2;
   [x,niter2]=ognewte(xi,x0,P0,A,BL,RA,CE,MITER,1/lam2,dtol,fast,flag);
   if flag>=1;disp('2nd point found');end
   xi2=x;
   niter=niter1+niter2;
   lam=lam2;tollam=tol*(lam1/lam2-1);dxn=10+tollam;indrun=0;
   [d2,dinf,Dtem]=hseldist(x,A,B,RA,CE,1,flag);
   xi_cor=0;
   
   while lam > lam0(2) & dxn(1) > tollam(1) & dinf>tol
      indrun=indrun+1;tloops=tloops+1;xi_cor=xi_cor-1;
%      d_fact=(dxn/(100*tollam/tol+dxn))^1.4;
      d_fact=(dinf/(1+dinf))^1.4;
%      d_fact=(d2/(1+d2))^1.6;
      lam=max(lam0(2),lam*lam_fact*d_fact);
      dinfi=0;xi_fact=1;
      while dinfi<=0
         atemp=(xi2-xi1)/(lam2-lam1);dellam=lam-lam2;
         xi=x+atemp*dellam*xi_fact;
         [d2,dinfi,Dtem]=hseldist(xi,A,B,RA,CE,1,flag);
         xi_fact=xi_fact/2;xi_cor=xi_cor+1;
      end
      BL=B+dinf*BLR*lam;
      [x,niter3]=ognewte(xi,x0,P0,A,BL,RA,CE,MITER,1/lam,dtol,fast,flag);
      [d2,dinf,Dtem]=hseldist(x,A,B,RA,CE,1,flag);
      niter=niter+niter3;
      xi1=xi2;xi2=x;lam1=lam2;lam2=lam;
      dxn=sqrt(abs((xi2-xi1)'*P0*(xi2-xi1)));
      tollam=tol*(lam1-lam2)/lam2;
      if flag>=1
         disp('     dinf      dx/tol_lam      n_iter        xi_cor         lam')
         disp([dinf dxn/tollam niter xi_cor lam])
      end
      
   end
   [x,niter3]=ognewte(x,x0,P0,A,BL,RA,CE,20,1/lam,dtol,fast,flag);
   niter=niter+niter3;
   nitert=[niter1,niter2,tloops,niter];
      if flag>=1
         disp('     dinf      dx/tol_lam      n_iter        xi_cor         lam')
         disp([dinf dxn/tollam niter xi_cor lam])
      end
end

