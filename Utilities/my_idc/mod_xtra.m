function [A,B,C,D,F,By,Dy,X0]=mod_xtr(FF,CC,thx,ninp,noutp,filord,nan_ind);

% function [A,B,C,D,F,By,Dy]=mod_xtr(f,q,thx,ninp,noutp,filord,nan_ind);
% extract and (A,B,C,D) representation from (f,q,theta)
%    also supplies the F,By matrices for the nonminimal
%    feedback representation (F,[B,By],C,[D,Dy]) (Dy=0)

[nth,ntem]=size(thx);nmax=max(filord);
if nargin<7;nan_ind=1;end
n_nan=length(nan_ind)-1;
n_thb=(ninp+2+n_nan)*nmax+ninp;

D=zeros(noutp,ninp);Dy=zeros(noutp,noutp);
A=[];B=[];C=CC;F=FF;By=[];X0=[];
   for i=1:noutp
      if i==1,nftem=0;else,nftem=sum(filord(1:i-1));end
      nftemi=nftem+filord(i);n=filord(i);
      By=[By,zeros(nftem,1);zeros(n,(i-1)),...
             thx(ninp*(1+n)+1:ninp+(ninp+1)*n,i)];
      B=[B;unvector(thx(1+ninp:ninp*(n+1),i),n,ninp)];
      D(i,:)=(thx(1:ninp,i))';
   end
B=B+By*D;A=F+By*C;

   for i_nan=1:n_nan+1
      xx0=[];
      for i=1:noutp
         n=filord(i);
         nan_loc=(ninp+1)*n+ninp+(i_nan-1)*n;
         xx0=[xx0;thx(nan_loc+1:nan_loc+n,i)];
      end
      X0=[X0,xx0];
   end

