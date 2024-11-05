function [xn,niter]=pe_newt(x0,tol,Nmax,cutoff,ff,q,cy,zq,RWi,N_Acns,t,u,y,ww);
%

xn=x0;errx=1;meth=2;arm=0;
stdes=0;indcall=1;
n=length(xn);df=0*xn;niter=0; 
ii=[1:n];iii=ii+n*(ii-1);
a=[1 .5 .1 .05 .01 1e-3 1e-4 1e-5 1e-6 1e-7 1e-8 1e-10 1e-12  0];
%a=[1 .5 .1 .05 .01 1e-3 1e-4 1e-5 1e-6 1e-7 1e-8 1e-10 1e-12 -1e-3 -1e-6 -1e-8 0];
la=length(a);  a_ind=1;

[fn,df,ddF]=pe_obj(xn,indcall,cutoff,ff,q,cy,zq,RWi,N_Acns,t,u,y,ww);

while errx > tol
   niter=niter+1;f=fn;
   if f>1.e50; disp('Optimization failure');return;end 
   x=xn;
   if stdes==0
         eif = svd(ddF);lamf=min(eif);lamff=max(eif);rcon=lamf/(1+lamff);
         if rcon(1) < 1.e-8
            ddF(iii)=ddF(iii)+(1.e-8*(lamff+1));
         end
         dk=-(ddF\df);
   else 
      zzz=abs(fn)/(1+df'*df); dk=-df*zzz;
   end
   % descent update
   a_ind=a_ind-1;if a_ind==0;a_ind=1;end
   xn=x+a(a_ind)*dk;indcall=0;

   [fn,df,ddF]=pe_obj(xn,indcall,cutoff,ff,q,cy,zq,RWi,N_Acns,t,u,y,ww);

   if fn>1.e50 | isnan(fn) ;err=1;else;err=fn-f;end
   while err>=0
      a_ind=a_ind+1;
      if a_ind>=la
         disp('armijo step at 0');err=-1;xn=x;fn=f;
      else
         xn=x+a(a_ind)*dk;  
         [fn,df,ddF]=pe_obj(xn,indcall,cutoff,ff,q,cy,zq,RWi,N_Acns,t,u,y,ww);
         if fn>1.e50 | isnan(fn) ;err=1;else;err=fn-f;end
      end
   end
   xn=x+a(a_ind)*dk;indcall=1;
   [fn,df,ddF]=pe_obj(xn,indcall,cutoff,ff,q,cy,zq,RWi,N_Acns,t,u,y,ww);
   errx=f-fn;
   %--------------------------------- Stopping criterion
     disp('O-Iteration,    Error,      Fun.value,    Stepsize,      x-xn')
     disp([niter errx fn a(a_ind) norm(x-xn)])
   if a(a_ind)==0
      if stdes==1
         errx=0;
         disp(': no cost reduction along descent direction');
      else
         stdes=1;errx=1;a_ind=1;
      end
   else
      stdes=0;a_ind=1;
   end
   if niter > Nmax; 
      errx=0;
      disp(': max no. of iterations exceeded');
   end
end


