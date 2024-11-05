function m=var_nan(x);

[N,M]=size(x);
m=zeros(1,M);
for i=1:M
  ind=find(~isnan(x(:,i)));
  m(i)=var(x(ind,i));
end
