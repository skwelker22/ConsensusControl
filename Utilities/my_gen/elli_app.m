function [c,R,Ri,vol]=elli_app(X,N);
clf
plot(X(1,:),X(2,:),'x');hold on

[n,meas]=size(X);
c=mean(X');c=c';
C=c*ones(1,meas);
re=max(sum(abs(C-X)));
R=eye(n,n)*re;Ri=eye(n,n)/re;
vol=zeros(N,1);k=0;

elliplot(Ri,c,3);pause(.1);
for i=1:N
    [i,N]
    k=k+1;
    d=rand(n,1)-.5;d=d/norm(d);elliplot(Ri,c,6);%plot(d(1),d(2),'gx');
%    d=X(:,i)-c;d=d/norm(d);elliplot(Ri,c,6);%plot(d(1),d(2),'gx');
    [c,R,Ri]=elli_cut(c,R,Ri,X,d);
    elliplot(Ri,c,3);pause(.1);
    vol(k)=det(R);
end

