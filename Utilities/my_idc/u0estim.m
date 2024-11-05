function [c0]=icestim(f,b,c,d,u,y,t,op);
%function c0=icestim(f,b,c,d,u,y,t,op);
% for the estimation of initial conditions of [f b c d]
% with i/o pairs [t,(u,y)].
% use all data for op=1

if nargin < 8;op=0;end
toler=1.e-6;
[nn,noutp]=size(y);[nn,nin]=size(u);nf=length(f);

    time_x=max(t);

ind_m=max(find(t<=time_x));txx=t(1:ind_m);

wu=lsim4(f,b,c,d,u(1:ind_m,:),txx);
Ferr=vector(y(1:ind_m,:)-wu);
BBB=zeros(ind_m*noutp,nin);
    ov=ones(size(txx));
    for i=1:nin
        ei=zeros(1,nin);ei(i)=1;
        w0=lsim4(f,b,c,d,ov*ei,txx);
        BBB(:,i)=vector(w0);
    end    
    c0=pinv([BBB])*Ferr;
