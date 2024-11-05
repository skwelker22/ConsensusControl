function x=vector(A);
% usage: x=vector(A)
% vectorize, in columns, a given matrix A

[n,m]=size(A);x=zeros(n*m,1);
for i=1:m
x(1+(i-1)*n:i*n)=A(1:n,i);
end
