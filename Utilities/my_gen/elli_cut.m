function [c,R,Ri]=elli_cut(c,R,Ri,X,d);
% function [c,R,Ri]=elli_cut(c,R,Ri,X,d);
% ellipsoid cuts to compute the min volume ellipsoid containing several vectors
% c,R,Ri: current estimate
% X: a matrix of the column vectors
% d: a direction (unit vector)

epsi=1e-4;

[n,meas]=size(X);
d=d/norm(d);
C=c*ones(1,meas);
prj=d'*(X);
pmin=min(prj);pmax=max(prj);
m=(pmin+pmax)*d/2; r=(pmax-pmin)/2;

ym=m'*d;
[c,R,Ri]=setmem(c,R,Ri,ym,d,r,epsi,1);

rmax=0;
for i=1:meas
    rr=(X(:,i)-c)'*Ri*(X(:,i)-c);rmax=max(rr,rmax);
end
Ri=Ri/rmax;R=R*rmax;vol=det(R);

for j=1:10
    dx=rand(n,1)-.5;dx=dx/norm(dx);
    c1=c+dx*norm(c)/100;rmax=0;
    for i=1:meas
        rr=(X(:,i)-c1)'*Ri*(X(:,i)-c1);rmax=max(rr,rmax);
    end
    Ri1=Ri/rmax;R1=R*rmax;vol1=det(R1);
    if vol1<vol; R=R1;Ri=Ri1;c=c1;vol=vol1;end
    
    dx=(rand(n,n)-.5);dx=dx/norm(dx);D=dx*dx'*norm(Ri)/10;
    Rix=Ri+D;
    for i=1:meas
        rr=(X(:,i)-c)'*Rix*(X(:,i)-c);rmax=max(rr,rmax);
    end
    Ri1=Rix/rmax;R1=inv(Ri1);R1=(R1+R1')/2;vol1=det(R1);
    if vol1<vol; R=R1;Ri=Ri1;vol=vol1;end
end