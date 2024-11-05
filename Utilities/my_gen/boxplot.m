function x=boxplot(A,B,xmin,xmax,ymin,ymax);

[m,n]=size(A);
x=0;sscol='bcgmrwy';colind=0;

x1=[xmin:(xmax-xmin)/100:xmax]';
y1=[ymin:(ymax-ymin)/100:ymax]';
axis([xmin,xmax,ymin,ymax])
for i=1:m;
colind=colind+1;
if colind >= 8;colind=1;end
hold on
if A(i,2)~=0
y=(B(i)-A(i,1)*x1)/A(i,2);
plot(x1,y,sscol(colind))
else
x=(B(i)-A(i,2)*y1)/A(i,1);
plot(x,y1,sscol(colind))
end
end