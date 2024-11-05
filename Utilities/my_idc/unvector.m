function A=unvector(x,n,m);
% usage: A=unvector(x,n,m)
% un-vectorize, in columns of length n, a vector x (m*n)

A=zeros(n,m);
for i = 1:m
A(:,i)=x((i-1)*n+1:i*n);
end
