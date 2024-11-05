function z=zoh(u,N);
% function z=zoh(u,N);
% zero order hold: holds the input value constant for N samples
% convolution with a window of 1's. Column vectors u

[nr,nc]=size(u);
h=ones(N,1);z=0*u;
for i=1:nc
    zt=conv(u(:,i),h);
    z(:,i)=zt(1:nr);
end

