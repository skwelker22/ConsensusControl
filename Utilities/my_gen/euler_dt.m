function Hd = euler_dt(H,T,fb)

if nargin <3,fb=0;end

H=ss(H);
if fb==0
    Hd=ss(H.a*T+eye(size(H.a)),H.b*T,H.c,H.d,T);
else
end