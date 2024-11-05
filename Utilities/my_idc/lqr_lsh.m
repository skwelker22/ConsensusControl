function [n_L,d_L]=lqr_lsh(num,den,lam,lp_s,tar);
% Usage:   [n_L,d_L]=lqr_lsh(num,den,lam,lp_s,tar);
% Compute lqr-based target for PID tuning

% K. Tsakalis 11/22/2001

if nargin<5;tar=0;end

w=logspace(log10(lam)-1,log10(lam)+3,4000);
%[A,B,C,D]=tf2ss(conv(num,[1,tar*lam]),conv([den 0],[1 100*lam]));
[A,B,C,D]=tf2ss(num,conv(den,[1 tar]));

Q=(C'*C)+lp_s(2)*A'*C'*C*A+lp_s(2)/100*sqrt(C*C')*eye(size(A));
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
R=sqrt(RM*RL); K=lqr(A,B,Q,R);  
[n_L,d_L]=ss2tf(A,B,K,D,1);
end

con=hi_ncf(ss(A,B,C,D));
loop=ss(A,B,C,D)*con;
[n_L,d_L]=tfdata(tf(loop),'v');