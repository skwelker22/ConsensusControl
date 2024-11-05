function [ap,bp,cp,dp]=pidss(kpid,T_der,T_fil,dt,pidform);
% [ap,bp,cp,dp]=pidss([kp,ki,kd],T_der,T_fil,dt,pidform);
% create a state-space representation of a PID
% if pidform = 1, then the summation form is assumed
% otherwise the pid gains are in the factored form.
% u=[kp+1/(ki*s)+kd*s/(T_der*s+1)][1/(T_fil*s+1)][e]
% dt=sampling time (0 for continuous time)

if nargin<5;pidform=0;end
if nargin<4;dt=0;end

kp=kpid(1);ki=kpid(2);kd=kpid(3);
[af,bf,cf,df]=tf2ss(1,[T_fil,1]);
if pidform==0
   [ai,bi,ci,di]=tf2ss(kp/ki,[1,0]);
   [ad,bd,cd,dd]=tf2ss([kp*kd,0],[T_der,1]);
elseif pidform==1
   [ai,bi,ci,di]=tf2ss(1/ki,[1,0]);
   [ad,bd,cd,dd]=tf2ss([kd,0],[T_der,1]);
elseif pidform==-1
   [ai,bi,ci,di]=tf2ss(ki,[1,0]);
   [ad,bd,cd,dd]=tf2ss([kd,0],[T_der,1]);
end

[ap,bp,cp,dp]=addss(ai,bi,ci,di,ad,bd,cd,dd);
dp=dp+kp;
if T_fil > 0
   [ap,bp,cp,dp]=series(af,bf,cf,df,ap,bp,cp,dp);
end
if dt >0
   [ap,bp,cp,dp]=bilin(ap,bp,cp,dp,1,'BwdRec',dt);
end
