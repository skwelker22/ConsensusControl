function rmsy=rms(Y);
%function rmsy=rms(Y);
% rms value of Y in columns, unless row vector

Y=abs(Y);
[n,m]=size(Y);
if n==1, Y=Y';end
[n,m]=size(Y);

rmsy=sqrt(sum(Y.*Y)/n);