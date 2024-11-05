function x=pid_eval(numo,deno,PID,FIL,lam,frw,fplantf);
% function x=pid_eval(numo,deno,PID,FIL,lam,frw,fplantf);
%   to evaluate the closed-loop pid performance/
% [numo,deno] describe the plant, PID contains the
%   PID numerator/denominator, FIL contains the prefilter,
%   lam is the desired BW and frw specifies the frequencies 
%   of interest

if nargin<7;fplantf=0;end

if fplantf==0
  x=0;
npid=PID(1,:);dpid=PID(2,:);nfil=FIL(1,:);dfil=FIL(2,:);
ztar=freqs(lam,[1 lam],frw);
zfil=freqs(nfil,dfil,frw);
clloop=conv(dpid,deno)+conv(npid,numo);
disp('closed loop poles')
rcl=roots(clloop);
disp(rcl)

sens=freqs(conv(dpid,deno),clloop,frw);
csen=freqs(conv(npid,numo),clloop,frw);
cont=freqs(conv(npid,deno),clloop,frw);
udis=freqs(conv(numo,dpid),clloop,frw);

clf;subplot(121)
loglog(frw,abs(sens),frw,abs(csen))
title('Cl.Lp S and T')
subplot(122)
loglog(frw,abs(cont),frw,abs(udis))
title('Cl.Lp T_r->u and T_du->y')
pause
clf;subplot(111)
loglog(frw,abs(csen.*zfil),frw,abs(ztar),'--')
title('Reference to output t.f.')
pause

tmax=100/lam;tstep=.1/lam;
t=[0:tstep:tmax];
y=step(conv(conv(npid,numo),nfil),conv(clloop,dfil),t);
y2=step(conv(dpid,numo),clloop,t);
index=find(abs(y-1)>1.e-2);index=[1:max(index)];
clf,subplot(121),plot(t(index),y(index))
title('Step resp. r->y')
subplot(122),plot(t(index),y2(index))
title('Step resp. ud->y')
pause
clf

else

npid=PID(1,:);dpid=PID(2,:);nfil=FIL(1,:);dfil=FIL(2,:);
zpid=freqs(npid,dpid,frw);
lpgn=zpid.*deno;
maglp=abs(lpgn);
phalp=phase(lpgn)*180/3.14159;
csen=lpgn./(1+lpgn);
sens=1.0./(1+lpgn);
cont=zpid./(1+lpgn);
udis=deno./(1+lpgn);
clf;subplot(121)
loglog(frw,abs(sens),frw,abs(csen))
title('Cl.Lp S and T')
subplot(122)
loglog(frw,abs(cont),frw,abs(udis))
title('Cl.Lp T_r->u and T_du->y')
pause


clf;subplot(121)
loglog(frw,abs(maglp));grid
title('Loop Gain Magnitude')
subplot(122)
semilogx(frw,phalp);grid
title('loop gain phase')
pause

gcross=min(find(maglp<=1));
fcross=frw(gcross);
phasmarg=phalp(gcross)+180;
disp('GAIN CROSSOVER,    GAIN,   PHASE MARGIN')
disp([fcross,maglp(gcross),phasmarg])

pcross=max(find(phalp>=-180));
if length(pcross)~=1;pcross=length(frw);end
fcross=frw(pcross);
gainmarg=1/maglp(pcross);
disp('PHASE CROSSOVER,    GAIN MARGIN,   PHASE')
disp([fcross,gainmarg,phalp(pcross)])

end