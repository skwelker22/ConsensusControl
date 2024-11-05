function [A,B,C,D,F,By,Dy,X0]=mod_xtr(FF,CC,thx,ninp,noutp,filord,nan_ind);

% function [A,B,C,D,F,By,Dy]=mod_xtr(f,q,thx,ninp,noutp,filord,nan_ind);
% extract and (A,B,C,D) representation from (f,q,theta)
%    also supplies the F,By matrices for the nonminimal
%    feedback representation (F,[B,By],C,[D,Dy]) (Dy=0)

[nth,ntem]=size(thx);nmax=max(filord);
if nargin<7;nan_ind=1;end
n_nan=length(nan_ind)-1;
n_ths=(ninp+2+n_nan)*nmax;
n_thb=(ninp+2+n_nan)*nmax+ninp;

   if nth>n_ths,nc=0;else,nc=1;end
   if nc==0
      if nth ~= n_thb
         disp('Incorrectly formed thx; reset to str.proper model')
         disp('Resulting model may be poor')
         nc=1;
      end
   end

D=zeros(noutp,ninp);Dy=zeros(noutp,noutp);
A=[];B=[];C=CC;F=FF;By=[];X0=[];
   for i=1:noutp
      if i==1,nftem=0;else,nftem=sum(filord(1:i-1));end
      nftemi=nftem+filord(i);n=filord(i);
%      f=F(nftem+1:nftemi,nftem+1,nftemi)';
%      q=C(i,nftem+1:nftemi)';
      By=[By,zeros(nftem,1);zeros(n,(i-1)),...
             thx(ninp*n+1:(ninp+1)*n,i)];
      B=[B;unvector(thx(1:ninp*n,i),n,ninp)];
        if nc==0
          D(i,:)=(thx((ninp+1)*n+1:(ninp+1)*n+ninp,i))';
        end
   end
B=B+By*D;A=F+By*C;

   for i_nan=1:n_nan+1
      xx0=[];
      for i=1:noutp
         n=filord(i);
         nan_loc=(ninp+1)*n+ninp*(1-nc)+(i_nan-1)*n;
         xx0=[xx0;thx(nan_loc+1:nan_loc+n,i)];
      end
      X0=[X0,xx0];
   end

