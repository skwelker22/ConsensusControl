function [P,K]=Kpid2tf(x,a);

if nargin<2;a=1;end

cs=0;P=[];K=[];
if isobject(x), cs=1;end
[n,m]=size(x); if n==2;cs=1;x=tf(x(1,:),x(2,:));end

if cs==0
    K=x;
    P=K(1)+tf(K(2),[1 0])+tf([K(3) 0],[K(4) 1]);
    if a~=1;
        P=ss(P);P.a=P.a/a;P.b=P.b/a;P=tf(P);
    end    
end

if cs==1
    P=x;K=zeros(1,4);
    if a~=1;
        P=ss(P);P.a=P.a/a;P.b=P.b/a;P=tf(P);
    end    
    [n,d]=tfdata(P,'v');
    d=d/d(2);
    if length(d)==2;K(4)=0;else;K(4)=d(1);end
    K(2)=n(3);
    K(1)=n(2)-K(2)*K(4);
    K(3)=n(1)-K(1)*K(4);
    tem=P;P=K;K=tem;
end
    