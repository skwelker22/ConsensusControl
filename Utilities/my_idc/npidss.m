function [a,b,c,d]=npidss(kpid,T_der,T_fil,dt,pidform);
% function [a,b,c,d]=npidss(kpid,T_der,T_fil,dt,pidform);
% create state-space representation of mimo PIDs (diagonal); 
% The i-th row of kpid should contain the gains of the
% i-th PID (kp,ki,kd);
% if pidform = 1, then the summation form is assumed
% u=[kp+1/(ki*s)+kd*s/(T_der*s+1)][1/(T_fil*s+1)][e]
% otherwise the pid gains are in the factored form.
% T_der, T_fil should be vectors with the corresponding TC's
% dt=sampling time (0 or [] for cont.time)

if nargin<5;pidform=0;end
if nargin<4;dt=0;end
if nargin<3;T_fil=0;end
[n,nx]=size(kpid);
if length(T_fil)<n;T_fil=T_fil(1)*ones(n,1);end
if length(T_der)<n;T_der=T_der(1)*ones(n,1);end
[a,b,c,d]=pidss(kpid(1,:),T_der(1),T_fil(1),dt,pidform);

if n>1
   for i=2:n
      [f,g,h,e]=pidss(kpid(i,:),T_der(i),T_fil(i),dt,pidform);
      [a,b,c,d]=append(a,b,c,d,f,g,h,e);
    end
end
