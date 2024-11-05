function [x0,c0]=icestim(f,b,c,d,u,y,t,op);
%function x0=icestim(f,b,c,d,u,y,t,op);
% for the estimation of initial conditions of [f b c d]
% with i/o pairs [t,(u,y)].
% use all data for op=1

if nargin < 8;op=0;end
toler=1.e-6;
[nn,noutp]=size(y);[nn,nin]=size(u);nf=length(f);

if op == 0
    lf_max=max(real(eig(f)));
    lama=-lf_max*0.85;
    PP=lyap(f+lama*eye(nf,nf),eye(nf,nf));
    sspp=eig(PP);l_ma=max(real(sspp));l_mi=min(real(sspp));
    x1=log(l_mi/l_ma)/2;x2=log(1+nf*trace(c*c'))/2;x3=log(toler);
    time_x=t(1)-1/(lama+1/2/l_ma)*(x1-x2+x3);
else
    time_x=max(t);
end

ind_m=max(find(t<=time_x));txx=t(1:ind_m);

wu=lsim4(f,b,c,d,u(1:ind_m,:),txx);
AAA=zeros(ind_m*noutp,nf);Ferr=vector(y(1:ind_m,:)-wu);
for i=1:noutp
    w0=lsim4(f',zeros(nf,1),eye(nf),zeros(nf,1),0*txx,txx,c(i,:)');
    AAA(1+(i-1)*ind_m:i*ind_m,:)=w0;
end
if nargout>1
    BBB=zeros(ind_m*noutp,nin);
    ov=ones(size(txx));
    for i=1:nin
        ei=zeros(1,nin);ei(i)=1;
        w0=lsim4(f,b,c,d,ov*ei,txx);
        BBB(:,i)=vector(w0);
    end    
    x0c=pinv([AAA BBB])*Ferr;
    x0=x0c(1:nf);c0=x0c(nf+1:nf+nin);
else
    x0=pinv(AAA)*Ferr;c0=[];
end
