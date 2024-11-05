function lab=labelme(t,y,t0);
% usage: lab=labelme(t,y,t0);
% labels multiple lines on the plot (t,y) at t=t0
% must be the current plot
% y is assumed to be plotted in columns
% t0 could be defined as a vector

lab=[];
[n,m]=size(y);
tx=zeros(1,m);
for i=1:m
	tx(i)=t0(min(i,length(t0)));
end
tin=tx;
for i=1:m
	ind=find (abs(t-tx(i))==min(abs(t-tx(i))));ind=ind(1);
	tin(i)=ind;
end
for i=1:m
text(t(tin(i)),y(tin(i),i),num2str(i))
end
