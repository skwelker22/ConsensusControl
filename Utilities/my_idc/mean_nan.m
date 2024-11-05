function m=mean_nan(x);

[N,M]=size(x);
m=zeros(1,M);
for i=1:M
  ind=find(~isnan(x(:,i)));
  m(i)=mean(x(ind,i));
end
