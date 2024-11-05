function [xn,niter]=ognewt(x0,xp,P0,A,B,RA,CE,Nmax,lam,tol,fast,flag);
%function [xn,niter]=ognewt(x0,xp,P0,A,B,RA,CE,Nmax,lam,tol,fast,flag);

xn=x0;errx=1;meth=2;arm=0;
stdes=0;indcall=1;
n=length(xn);df=0*xn;niter=0; 
ii=[1:n];iii=ii+n*(ii-1);
a=[1 .5 .1 .05 .01 1e-3 1e-4 1e-5 1e-6 1e-7 1e-8 1e-9 1e-10 1e-12 1e-14 1e-16 0];
la=17;  a_ind=1;

[fn,df,ddF]=proj_obj(xn,xp,P0,A,B,RA,CE,lam,indcall);

while errx(1) > tol(1)
   niter=niter+1;f=fn;
   if f(1)>1.e50;if flag(1)>1; disp('Optimization failure');end;return;end 
   x=xn;
   if stdes==0
      if fast(1)>0
         eif = svd(ddF);lamf=min(eif);lamff=max(eif);
         rcon=lamf/(1+lamff);
         if rcon(1) < tol/10
            ddF(iii)=ddF(iii)+(tol/10*(lamff+1));
            if flag(1)>1;
               disp(['ognewt: hessian correction: cond.num. = ',num2str(rcon)]);
            end
         end
         dk=-(ddF\df);
      elseif fast(1)==0
         ddF(iii)=ddF(iii)+tol/10;
         dk=-(ddF\df);
      else
         dk=-(ddF\df);
      end
   else 
      zzz=abs(fn(1))/(1+df'*df);
      dk=-df*zzz;
   end
   % descent update
   %   NNd=sqrt(abs(df'*dk));
   %   if NNd(1)<=0.25, a=a0;else,a=a0/(a0+NNd);end
   a_ind=a_ind-1;if a_ind==0;a_ind=1;end
   xn=x+a(a_ind)*dk;indcall=0;
   [fn,df,ddF]=proj_obj(xn,xp,P0,A,B,RA,CE,lam,indcall);
   if fn(1)>1.e50;err=1;else;err=fn-f;end
   while err(1)>=0
      a_ind=a_ind+1;
      if a_ind>=la
         if flag(1)>1;disp('armijo step at 0');end
         err=-1;xn=x;fn=f;
      else
         xn=x+a(a_ind)*dk;  
         [fn,df,ddF]=proj_obj(xn,xp,P0,A,B,RA,CE,lam,indcall);
         if fn(1)>1.e50;err=1;else;err=fn-f;end
      end
   end
   xn=x+a(a_ind)*dk;indcall=1;
   [fn,df,ddF]=proj_obj(xn,xp,P0,A,B,RA,CE,lam,indcall);
   errx=f-fn;
   %--------------------------------- Stopping criterion
   if flag(1) > 3
      disp('O-Iteration,    Error,      Fun.value,    Stepsize,      Descent')
      disp([niter errx fn a(a_ind) stdes])
   end
   if a(a_ind)==0
      if stdes==1
         errx=0;
         if flag(1)>1;disp('ognewt: no cost reduction along descent direction');end
      else
         stdes=1;errx=1;a_ind=1;
      end
   else
      stdes=0;a_ind=1;
   end
   if niter > Nmax(1); 
      errx=0;
      if flag(1)>1;disp('ognewt: max no. of iterations exceeded');end
   end
end

