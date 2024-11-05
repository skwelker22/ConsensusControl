function [n_L,d_L]=lqr_lshs(A,B,C,D,lam,lp_s,tar);
% Usage:   [n_L,d_L]=lqr_lshs(A,B,C,D,lam,lp_s);
% Compute lqr-based target for PID tuning

if nargin<5;tar=1;end

w=logspace(log10(lam)-1,log10(lam)+3,4000);
%[A,B,C,D]=tf2ss(conv(num,[1,tar*lam]),conv([den 0],[1 100*lam]));
%[A,B,C,D]=tf2ss(num,[den 0]);
[ax,bx,cx,dx]=tf2ss(1,[1 0]);
[A,B,C,D]=series(ax,bx,cx,dx,A,B,C,D);
Q=C'*C+lp_s(2)*sqrt(C*C')*eye(size(A));
N=C'*D;D2=D'*D;
if lp_s(3) > 0; 
  R=lp_s(3);K=lqr(A,B,Q,R);
  [n_L,d_L]=ss2tf(A,B,K,D,1);
else
  R=1;
     K=lqr(A,B,Q,R);[n1,d1]=ss2tf(A-B*K,B,K,-1,1);m1=bode(n1,d1,w);
     ind=max(find(m1<0.707));if isempty(ind);ind=1;end;bw=w(ind);
     while bw>lam
       R=R*100;
       K=lqr(A,B,Q,R);[n1,d1]=ss2tf(A-B*K,B,K,-1,1);m1=bode(n1,d1,w);
       ind=max(find(m1<0.707));if isempty(ind);ind=1;end;bw=w(ind);
     end
  RM=R;
     K=lqr(A,B,Q,R);[n1,d1]=ss2tf(A-B*K,B,K,-1,1);m1=bode(n1,d1,w);
     ind=max(find(m1<0.707));if isempty(ind);ind=1;end;bw=w(ind);
     while bw<lam
       R=R/100;
       K=lqr(A,B,Q,R);[n1,d1]=ss2tf(A-B*K,B,K,-1,1);m1=bode(n1,d1,w);
       ind=max(find(m1<0.707));if isempty(ind);ind=1;end;bw=w(ind);
     end
  RL=R;
     ER=log(RM)-log(RL);
  while abs(ER)>0.01 
     R=sqrt(RM*RL);
     K=lqr(A,B,Q,R);[n1,d1]=ss2tf(A-B*K,B,K,-1,1);m1=bode(n1,d1,w);
     ind=max(find(m1<0.707));if isempty(ind);ind=1;end;bw=w(ind);
     if bw>lam;RL=R;else;RM=R;end
     ER=log(RM)-log(RL);
  end
R=sqrt(RM*RL); K=lqr(A,B,Q,R);  [n_L,d_L]=ss2tf(A,B,K,D,1);
end
