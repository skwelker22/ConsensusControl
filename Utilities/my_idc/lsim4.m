function  yout = lsim4(a, b, c, d, u, t, x0)
t=[0:length(t)-1]'*(t(2)-t(1));
if nargin<7;
yout=lsim(a,b,c,d,u,t);
else
yout=lsim(a,b,c,d,u,t,x0);
end

