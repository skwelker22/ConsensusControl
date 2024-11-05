function x=elliplot(RA,CE,colind);
% function x=elliplot(RA,CE,colind);
% plots the ellipsoid (x-CE)'RA(x-CE) < 1

[n,m]=size(CE);

if n==2

x=0;sscol='bcgmrwy';
if nargin<3;colind=0;end

w=[0:5:360]*3.14159/180;
unitvec=[sin(w);cos(w);zeros(n-2,1)];

for i=1:m;
colind=colind+1;
if colind >= 8;colind=1;end
cent=CE(:,i);
radi=RA(:,(i-1)*n+1:i*n);
rri=inv(sqrtm(radi));
xx=rri*unitvec;
hold on
plot(xx(1,:)+cent(1),xx(2,:)+cent(2),sscol(colind))
plot(cent(1),cent(2),[sscol(colind),'+'])
end

end
    