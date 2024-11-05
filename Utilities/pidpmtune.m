function [pid,cpid,dpid]=pidpmtune(bw,g,tau,pm,n,T, pole)

% [pid,cpid,dpid]=pidpmtune(bw,g,tau,pm,n,T)

% bw = bandwidth

% g = system

% tau = derivative TC

% pm = phase margin

% n = PI/PID order (force)

% T = sample time

% pole = -rhp pole

 

gc=bw/1.5;


if nargin<3,tau=.01; end

if isempty(tau),tau=0.1/gc,end

if nargin<4;pm=50;end

if isempty(pm);pm=50;end

if nargin<5;n=[];end

if nargin<6;T=0;end

if nargin<7;pole=0;end

dpid=[];


if pole == 0

cpid=tf(1,[1,0]);

[m,p]=bode(g*cpid,gc);

p=mod(p,360);

if p>0;p=p-360;end

th=(-180-p+pm)

if isempty(n)

    if abs(th)>75; n=2;else;n=1;end

end

else

    n=2;

end


if n==1;cpid=tf(1,[1 0]);else;cpid=tf(1,conv([tau,1],[1 -pole]));end

[m,p]=bode(g*cpid,gc);

p=mod(p,360);

if p>0;p=p-360;end

th=(-180-p+pm)

tz=tan(max(th,0)/n*pi/180)/gc

if n==2;cpid=tf([tz*tz 2*tz 1],conv([tau,1],[1 -pole]));else;cpid=tf([tz 1],[1,0]);end

[m,p]=bode(g*cpid,gc)

k=1/m

cpid=cpid*k;


[nu,de]=tfdata(cpid,'v');


if length(nu)<3;nu=[0,nu];end

pid=[nu(2)-nu(3)*tau,nu(3),nu(1)-tau*(nu(2)-nu(3)*tau)];

if T>0,dpid=c2d(cpid,T,'tustin');end