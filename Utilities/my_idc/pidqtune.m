function [PID,Ub,Cpid]=pidqtune(lam,num,den,lp_s,f_,fr_r,flag);
%  Usage:
%  [PID,Ub]=PIDqtune(lam,numo,deno,[lp_s,f_,fr_r,flag]);
%
%  Computes the PID gains kp,ki,kd, by loop-shaping techniques to achieve
%     a closed-loop sensitivity bandwidth approx. equal to lam,
%     while maintaining a small sensitivity peak value.
%  [Ref. Grassi & Tsakalis, IEEE Trans. Contr. Sys. Techn., Sep. 2000, pp 842-847]
%
%  PID in additive form: u = [kp+ki/(s+lp_s(4))+kd*s/[T*s+1]] e
%     lp_s(4) can modify the integrator to a RHP pole (useful in cart-pendulum)
%
% Inputs:
%   lam: desired bandwidth; constrained by robustness spec's.
%   num,den: numerator,denominator of nominal plant.
%     IN V.6: lam can be a system object target loop and num a system
%     object plant; den = []; 
%   lp_s:  Row vector specifying the target loop shape via an LQR solution.
%      lp_s(1)= 1+r: unused here ---- (L2 optimization with Linf
%      constraint)
%      lp_s(2)= integral action penalty (can be used to adjust overshoot; def. 1e-3)
%      lp_s(3)= control input penalty (if zero, iterate until S-bandwidth=lam)
%      lp_s(4)= "integrator" corner frequency (-ve for controller with RHP pole)
%   f_(1): filter time constant for the pseudo-derivative term; default = 0.01;
%   f_(2): 0 for PID, 1 for PI
%   fr_r: Frequency range for the optimization;
%         =0: 300 pts in 0.01*lam--100*lam;
%         or, it may contain a row vector of user specified frequencies.
%         if so, target loop can be specified in fr_r(2,:).
%
%  Zero entries select default values; incomplete vectors are padded with zeros.
%
%  Intermediate Outputs (controlled by flag):
%  - Ub: fitting error; must be less than 1 to guarantee closed loop
%          stability; typical values should be around 0.2-0.3,
%  - PID zeros and Closed-loop poles
%  Displays:
%       Closed loop sensitivity, complementary sensitivity,
%       control sensitivity and response to step reference
%       commands and disturbances
%  Display suppressed when flag = -1
%
% OUTPUT PID format: [ numerator polynomial;
%                      denominator polynomial;
%                      kp, ki, kd;]
%
%  Limitations: The plant must be stabilizable by a PID with positive gains.
%     Scaling may be required if the optimizer magnitude exceeds the
%     radius of the default search set.
%

% K. TSAKALIS, 8/23/01 (quick version)
%      rev.    Oct 2004

TL=[];
temp=version;
%if str2num(temp(1))>=6;
    if isobject(num);PLANT=num;[num,den]=tfdata(tf(num),'v');end
    if isobject(lam);
        w=logspace(-8,8,1000)';
        TL=lam;
        [m]=bode(TL,w);m=m(:);figure,loglog(w,m);
        item=max(find(m>=0.707));lam=w(item);
        w=logspace(log10(lam/10),log10(lam*10),1000);
        [m]=bode(TL,w);m=m(:);figure,loglog(w,m);
        item=max(find(m>=0.707));lam=w(item);
    end
%end

FIL=[1;1];out_flag=0; con_=1e8; toler=1e-6;
format short e
while length(num) < length(den),num=[0 num];end

if nargin<4, lp_s=0; end
if nargin<5, f_=0;end
if nargin<6, fr_r=0;end
if nargin<7, flag=1;end

while length(lp_s)<4;lp_s=[lp_s,0];end
if f_(1)==0;f_(1)=0.01;end;T=f_(1);
if length(f_)<2;f_=[f_,0];end; c_=f_(2);
if lp_s(2)==0;lp_s(2)=0.001;end
op_t=lp_s(1);

[Nfrr,Mfrr]=size(fr_r);
if length(fr_r)==1
    fmin=0.01*lam;fmax=100*lam;
    frw=logspace(log10(fmin),log10(fmax),300);
end
if length(fr_r)>2;frw=fr_r(1,:);end
if isempty(TL)
    if Nfrr==1
        tar=lp_s(4);
        [n_L,d_L]=lqr_lsh(num,den,lam,lp_s,tar);
        lshape=freqs(n_L,d_L,frw).';
    else
        lshape=fr_r(2,:);
    end
else
    [n_L,d_L]=tfdata(TL,'v');
    lshape=freqs(n_L,d_L,frw).';
end
% create frequency responses
pla=freqs(num,den,frw).';
if flag >= 1
    %clf
    figure,
    loglog(frw,abs(lshape),frw,abs(pla));
    legend('loop','plant');
    title('Target loop and plant');pause(1)
end

%----- DEFINE FILTERS
dpid=conv([T,1],[1,lp_s(4)]);
npid1=[0 0 1];npid2=[0 1 0];npid3=[1 0 0];
zp1=freqs(npid1,dpid,frw).';
zp2=freqs(npid2,dpid,frw).';
zp3=freqs(npid3,dpid,frw).';
wei=1./(1+lshape);zi=lshape.*wei;
w1=zp1.*pla.*wei;w2=zp2.*pla.*wei;w3=zp3.*pla.*wei;
W=[w1,w2,w3];
%------ Linf - L2 Procedure
%       op_t = 0
[K,Ub]=lmiqpid(W,zi,toler,con_,T,c_,0,frw',flag);K
if c_==1;K(3)=T*K(2)-T*T*K(1);end
WK=W*K-zi;


if norm(K)>=con_(1)
    out_flag=1;
    if flag>=0
        disp('PIDQTUNE WARNING 1: Optimizer outside initial search set (con_)')
        disp('Optimizer magnitude and search set radius')
        disp([norm(K),con_(1)])
    end
end
if flag >=0
    %clf;
    figure,loglog(frw(1,:),abs(WK));
    title('Approximation Error');pause(1)
    disp(['Approximation error ',num2str(Ub)])
end

if op_t > 1
    op_t=max(op_t,1.001);
    op_ty=-op_t*Ub;
    [K,Ub]=lmiqpid(W,zi,toler,con_,T,c_,op_ty,frw',flag);
    if c_==1;K(3)=T*K(2)-T*T*K(1);end
    WK=W*K-zi;
    if norm(K)>=con_(1)
        out_flag=1;
        if flag>=0
            disp('PIDQTUNE WARNING 2: Optimizer outside initial search set (con_)')
            disp('Optimizer magnitude and search set radius')
            disp([norm(K),con_(1)])
        end
    end
    if flag >=0
        %clf;
        figure,loglog(frw(1,:),abs(WK));
        title('Approximation Error');pause(1)
        disp(['Approximation error ',num2str(Ub)])
    end
end

% Compute the pid gains
npid=K(1)*npid1+K(2)*npid2+K(3)*npid3;
if lp_s(4)==0
    ki=K(1);kp=K(2)-T*ki;kd=K(3)-kp*T;
else
    [R,P,kpp]=residue(npid,dpid);
    ki=R(2);    kd=-R(1)/P(1)/P(1);   kp=kpp-R(1)/P(1);
end
if flag >=2
    disp('PID zeros')
    disp(roots(npid))
end
while length(npid) < length(dpid),npid=[0 npid];end
PID=[npid;dpid;kp,ki,kd];

if flag >=1
    x=pidfeval(num,den,PID,0,FIL,lam,frw,0,flag);
end
Cpid=tf(npid,dpid);

%--------------------------------------------------------------------------

function [K,Ub]=lmiqpid(W,zi,toler,con_,T,c_,op_t,fr,flag);
% function [K,Ub]=LMIqPID(W,zi,toler,con,T,c_,op_t,fr,flag);
%  Convex optimization function to perform PID tuning

% K. TSAKALIS, 8/23/01

flag =0;
K=[0;0;0]; Aell=con_(1)*con_(1)*eye(length(K),length(K));
L=[];U=[];err=1;erro=1;kount=0;phi=0;kkount=0;indp=0;
op_tx=op_t;if op_tx(1) < 0;op_tx(1)=1;end

while err>toler
    [h,psi]=pidqcns(K,toler/10,T,c_);
    if norm(h) == 0 & op_t(1) < 0
        [h,psi]=pidqobj(W,zi,K,op_t,fr);
    end
    h0=norm(h);
    if h0>0
        % 'constraint-iteration'
        iterty=-1;deepcut=psi;
    else
        % 'objective-iteration'
        iterty=1;
        [h,phi]=pidqobj(W,zi,K,op_tx,fr);
        erro=sqrt(h'*Aell*h);deepcut=0;
        %     L=[L;phi-erro];U=[U;phi];
        if length(U)<1;U=phi;else;U=min(U,phi);end
    end
    [K,Aell,err_flag]=lmi_upd(K,h,Aell,deepcut,flag);
    kount=kount+1;kkount=kkount+1;
    h0=norm(pidqcns(K,toler/10,T,c_));
    if op_t(1) < 0
        [hx,h1]=pidqobj(W,zi,K,op_t,fr);
    else
        h1=0;
    end
    err=(h0)+h1+erro;
    if err_flag<0
        if flag>=0;disp('LMIqPID WARNING: Optimization failure');end
        err=0;
    end
end
if err_flag<0
    Ub=inf;
else
    Ub=sqrt(U);
end
%--------------------------------------------------------------------------

function  [h,val]=pidqcns(K,tol,T,c_);
%  function  [h,val]=PIDqcns(K,tol,T,c_);
%    specifying the PID parameter constraints
%

% K. TSAKALIS, 8/11/04

cc=inf;
if c_==1;cc=2*tol;end
h=0*K; val=0;
if K(1)<=tol,h(1)=-1;val=tol-K(1);
elseif K(2)<=tol,h(2)=-1;val=tol-K(2);
elseif K(3)<=tol,h(3)=-1;val=tol-K(3);
end
if K(3)-T*K(2)+T*T*K(1) >= cc
    h=[T*T, -T, 1]';
    val=(K(3)-T*K(2)+T*T*K(1))-cc;
end
%--------------------------------------------------------------------------

function [h,phi]=pidqobj(W,zi,K,op_t,fr);
% function [h,phi]=PIDqOBJ(W,zi,K,op_t,fr);
% Defines the optimization objective in pidftune

% K. TSAKALIS, 8/11/04

temp = W*K-zi;npt=length(temp);
if op_t(1)==1
    dw   = [fr;0]-[0;fr];dw=dw(1:npt);
    ph   = (abs(temp)).^2;
    phi  = (ph(1:npt-1)+ph(2:npt))'*dw(2:npt)/2+ph(1)*dw(1);
    h    = (temp(1:npt-1).*dw(2:npt))'*W(1:npt-1,:);
    h    = h+(temp(2:npt).*dw(2:npt))'*W(2:npt,:);
    h    = real(h+2*temp(1)'*W(1,:)*dw(1));
    h    = h'/norm(h);
elseif op_t(1) < 0
    phim = max(abs(temp))^2;
    phi  = max(abs(temp))^2 - op_t(1)*op_t(1);
    if phi > 0
        indp = min(find(abs(temp)==sqrt(phim)));
        h    = 2*(real(temp(indp)'*W(indp,:)))';
    else
        phi = 0; h = 0*K;
    end
else
    phi  = max(abs(temp))^2;
    indp = min(find(abs(temp)==sqrt(phi)));
    h    = 2*real(temp(indp)'*W(indp,:));
    h    = h'/norm(h);
end
%--------------------------------------------------------------------------

function [ncent,nPcov,err_flag]=lmi_upd(cent,subgr,Pcov,deep_cut,flag);
%  Usage: [ncent,nPcov,err_flag]=lmi_upd(cent,subgr,Pcov,deep_cut,flag);

% K. TSAKALIS, 8/23/96

zero_tol=1.e-8;err_flag=0;
ncent=cent;nPcov=Pcov;n=length(cent);
Inmat=eye(n,n);
if  nargin < 4,deep_cut=0;end
if  nargin < 5,flag=0;end
if deep_cut < 0
    if flag>=0
        disp('ERROR (LMI_UPD):')
        disp('Invalid Constraint Spec')
    end
    err_flag=-3;
    return
end
s=(subgr'*Pcov*subgr);
if s<-zero_tol
    if flag>=0
        disp('ERROR (LMI_UPD):')
        disp('Invalid Ellipsoid Matrix')
    end
    err_flag=-1;return;
elseif s<=0
    if flag>=0
        disp('WARNING (LMI_UPD):')
        disp('Ellipsoid thickness is zero along gradient'),
    end
    Pcov=Pcov+Inmat*zero_tol;
    s=(subgr'*Pcov*subgr);
    deep_cut=0;err_flag=1;
end
s=sqrt(s);g=subgr/s;a=deep_cut/s;
if a>1,
    if flag>=0
        disp('ERROR (LMI_UPD): Empty Feasible Set')
    end
    err_flag=-2;
    return
end
ncent=cent-(1+n*a)*Pcov*g/(n+1);
nPcov=(Pcov-2*(1+n*a)/(n+1)/(1+a)*Pcov*g*g'*Pcov)*n*n/(n*n-1)*(1-a*a);
%--------------------------------------------------------------------------

function [n_L,d_L]=lqr_lsh(num,den,lam,lp_s,tar);
% Usage:   [n_L,d_L]=lqr_lsh(num,den,lam,lp_s,tar);
% Compute lqr-based target for PID tuning

% K. Tsakalis 11/22/2001

if nargin<5;tar=0;end

w=logspace(log10(lam)-1,log10(lam)+3,4000);
%[A,B,C,D]=tf2ss(conv(num,[1,tar*lam]),conv([den 0],[1 100*lam]));
[A,B,C,D]=tf2ss(num,conv(den,[1 tar]));

Q=(C'*C)+lp_s(2)*A'*C'*C*A+lp_s(2)/100*(C*C')*eye(size(A));
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


%--------------------------------------------------------------------------

function [x,G]=pidfeval(numo,deno,PID,TF,FIL,lam,frw,fplantf,flag);
% function [x,G]=pidfeval(numo,deno,PID,TF,FIL,lam,frw,fplantf,flag);
%   to evaluate the closed-loop pid performance/
% [numo,deno] describe the plant, PID contains the
%   PID numerator/denominator, FIL contains the prefilter,
%   lam is the desired BW and frw specifies the frequencies
%   of interest

% K. TSAKALIS, 8/23/96

if nargin<8;fplantf=0;end
if nargin<9;flag=0;end
x=0;G=0;
npid=PID(1,:);dpid=PID(2,:);nfil=FIL(1,:);dfil=FIL(2,:);
if TF>1.e-4,npid=[0,npid];dpid=conv(dpid,[TF,1]);end
ztar=freqs(lam,[1 lam],frw);
zfil=freqs(nfil,dfil,frw);
clloop=conv(dpid,deno)+conv(npid,numo);
if flag>=1
    disp('closed loop poles')
    rcl=roots(clloop);
    disp(rcl)
end
sens=freqs(conv(dpid,deno),clloop,frw);
csen=freqs(conv(npid,numo),clloop,frw);
cont=freqs(conv(npid,deno),clloop,frw);
udis=freqs(conv(numo,dpid),clloop,frw);
G=freqs(conv(npid,numo),conv(dpid,deno),frw);
%clf;
if flag>=1
    figure,
    loglog(frw,abs(sens),frw,abs(csen));grid
    title('Cl.Lp S and T')
    pause(1)
end
if flag>=2
    figure,
    loglog(frw,abs(cont),frw,abs(udis));grid
    title('Cl.Lp T_r->u and T_du->y')
    pause(1)
end
if flag >=3
    %clf;
    figure,
    loglog(frw,abs(csen.*zfil),frw,abs(ztar),'--');grid
    title('Reference to output t.f.')
    pause(1)
end
tmax=15/lam;tstep=.1/lam;
t=[0:tstep:tmax];
y=step(conv(conv(npid,numo),nfil),conv(clloop,dfil),t);
y2=step(conv(dpid,numo),clloop,t);
if flag>=1
    %clf,
    figure,subplot(121),plot(t,y);grid
    title('Step resp. r->y')
    subplot(122),plot(t,y2);grid
    title('Step resp. ud->y')
    pause(1)
end
x=[t',y,y2];
