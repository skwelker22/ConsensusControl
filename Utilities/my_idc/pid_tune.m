function [kp,ki,kd,PID,WK,lshape]=PID_tune(lam,numo,deno,T,B,frw,lw,pos);
%Usage: [kp,ki,kd,PID]=PID_tune(lam,numo,deno,[T,B,frw,lw,pos]);
%
%  Computes the PID gains kp,ki,kd, tuned to
%  achieve a closed-loop bandwidth approx. equal to lam,
%  while maintaining a 'small' sensitivity peak value.
%
%  Limitations: The plant must be stabilizable by a PID with 
%      positive gains. Also, suitable choice of the target loop 
%      may be required for plants with integral action to achieve
%      adequate input disturbance attenuation.
%
% PID definition: u = kp*(1+1/[ki*s]+kd*s/[T*s+1]) e
%
% Inputs:  
%    lam: desired bandwidth; constrained by robustness
%        spec's and modelling errors; 
%    numo,deno: numerator,denominator of nominal plant.
%        The plant can also be specified by its frequency response;
%        in this case, set numo=0 and deno=row vector of freq.resp.
%        evaluated at the frequency points in frw (e.g., by freqs);
%        also, for this case, the option B(3)=3 must be used. 
%    T:  filter time constant for the pseudo-derivative term;
%        typically 0.5 or 1 * sampling time; default = 0.01
%    B:  vector specifying search set constraints and target loop:
%        B(1)=radius, 
%        B(2)=high-frequency gain constraint s.t. 
%          cl-lp input sensitivity < B(2)/T (default 1e3), 
%        B(3)=0 selects a single integrator as a target loop 
%            (OK for plants with no poles near the jw-axis)
%           ~=0 selects high order target loop with pole/zero 
%            specified by B(4,5) according to the following 
%            convention: 1-3: L(s)=lam(s+a)/(s+p), 4: L(jw) in lam
%        B(3)=1: a=lam*B(4), p=smallest in magnitude plant pole; 
%            here larger values for B(4) may induce sensitivity peaking;
%            smaller values may affect input disturbance attenuation; 
%            typical choices are around 0.2
%        B(3)=2: a, p are the smallest in magnitude plant pole and zero. 
%            (OK for plants with a pole-zero pair near the jw-axis)
%        B(3)=3: a=B(4);p=B(5)
%        B(3)=4: L(s) in frw(2,:); in this case lam is the expected BW 
%    default B=[1.e2,1.e3,0,0,0].
%    frw: row vector containing frequency points; default: 500
%          points equally spaced in log-scale between 
%          0.02*min|root(numo,deno)| 50*max|root(numo,deno)|
%        it can also be used to specify the frequency
%        response of the loop shape at the frequencies frw 
%    lw: frequency cut-off for the fitting weight (ratio to lam);
%        if <1, low frequencies are given higher weight; 
%        if >1, high frequencies are given higher weight; 
%        default: 0.1; 
%        the fitting weight can also be specified by using
%        lw(1,:)=num_weight;lw(2,:)=den_weight
%    pos: 0 restricts pid gains to be positive; 1 selects the
%          alternative PID kp+1/[ki*s]+kd*s/[T*s+1]
%          restricted to be minimum phase only)
%  zero entries select default values
%
%  Intermediate Outputs:
%  - Progress of Convex Optimization
%  - Ub: fitting error; must be less than 1 to guarantee
%         closed loop stability; typical value should be
%         around 0.2,otherwise large sensitivity peaks may occur; 
%         if too high, a suitable scaling can be achieved by
%         reducing lam and kp; note: kp~lam~Ub; tuning needs 
%         to be repeated only when lam should be scaled by
%         more than  an order of magnitude.)
%  - Resulting magnitude of optimizer; if larger than B(1), 
%         increase B(1) and recompute
%  - Closed-loop poles
%  
%  Displays:
%       Closed loop sensitivity, complementary sensitivity,
%       control sensitivity and response to step reference
%       commands and disturbances
%
%  Calls to LMI_PID_tune, PID_constraints, LMI_update;
%     PID parameter constraints may be specified through 
%     PID_constraints;  default constraints are that the PID
%     parameters are positive; in an alternative configuration, 
%     the continuous time PID transfer function must be minimum
%     phase, which does not exclude the possibility that some of
%     the PID parameters may be negative; in this case (pos=1)
%     the PID gains are returned in the alternative form.
%     Also, the initial search set radius is set by LMI_PID_tune
%     to B(1). Larger values may slow down convergence, 
%     while smaller values may overly restrict the PID parameters.
%     

FIL=[1;1];al=[];pm=[];
flag=1;     % set flag =0 to disable plots

if nargin<4, pos=0;lw=0;frw=0;B=0;T=0;
   elseif nargin<5, pos=0;lw=0;frw=0;B=0;
   elseif nargin<6, pos=0;lw=0;frw=0;
   elseif nargin<7, pos=0;lw=0;
   elseif nargin<8, pos=0;
end
if T == 0; T = 0.01;end
if B == 0; B=[0 0 0 0 0];end
if B(1) == 0;B(1)=1.e3;end
if B(2) == 0;B(2)=1.e3;end
if lw == 0; lw = 0.1;end
if frw == 0
   r1=abs(roots(numo));r2=abs(roots(deno));
   fmin=([r1;r2]);fmin=fmin(find(fmin>0));
   fmin=0.02*min(fmin);fmax=50*max([r1;r2]);
   fmax=max(10/T,fmax); fmin=min(fmin,.01*lam);
   frw=logspace(log10(fmin),log10(fmax),500);
end
[nnls1,nnls2]=size(frw);
if nnls1 >1 B(3)=4;end
if numo==0; fplantf=1; else, fplantf=0;end

if fplantf ==0
r2=roots(deno);pmin=r2(find(abs(r2)==min(abs(r2))));pmin=real(pmin(1));
r2=roots(numo);
if length(r2)>0;zmin=r2(find(abs(r2)==min(abs(r2))));zmin=real(zmin(1));
else, zmin=[];end
end
if B(3)==1,al=B(4)*lam;pm=-pmin;end
if B(3)==2,al=-zmin;pm=-pmin;end
if B(3)==3,al=B(4);pm=B(5);end
BDZ=0;

% create frequency responses
  dpid=[T 1 0];npid1=[0 1 0];npid2=[0 0 1];npid3=[1 0 0];
  zp1=freqs(npid1,dpid,frw(1,:)).';zp2=freqs(npid2,dpid,frw(1,:)).';
  zp3=freqs(npid3,dpid,frw(1,:)).'; 
     if fplantf==0
       pla=freqs(numo,deno,frw(1,:)).';
     else
       pla=deno.';
     end
     if B(3) == 0
       lshape=freqs(lam,[1 0],frw(1,:)).'; 
     elseif B(3)~= 4
       lshape=freqs(lam*[1 al],[1 pm 0],frw(1,:)).'; 
          if flag==1
            [lshapen,lam,al,pm]=lpshptun(frw(1,:),lshape,lam,al,pm);
            lshape=lshapen;
          end
     else
       lshape=(frw(2,:)).';
     end
disp('Loop shape parameters lam -zero -pole')
disp([lam,al,pm])
[nlw,mlw]=size(lw);
     if length(lw)==1
       if lw<1
        lwwei=(freqs([1 lam],[1 lw*lam],frw(1,:))).';
       else
        lwwei=(freqs([1/lam 1],[1/(lam*lw) 1],frw(1,:))).';
       end
     elseif min(nlw,mlw)==2
       lwwei=(freqs(lw(1,:),lw(2,:),frw(1,:))).';
     else 
       lwwei=lw;
     end
  wei=lwwei./(1+lshape);zi=wei.*lshape;        w1=zp1.*pla.*wei;w2=zp2.*pla.*wei;w3=zp3.*pla.*wei;
  W=[w1,w2,w3];
  wei0=(1)./(1+lshape);zi0=wei0.*lshape;      
  W0=[zp1.*pla.*wei0,zp2.*pla.*wei0,zp3.*pla.*wei0];

% Use convex optimization
  [K,Ub]=LMI_PID(W,zi,5.e-5,T,pos,B,BDZ);
WK=W0*K-zi0;
disp('Approximation error')
disp(Ub)
disp('Optimizer magnitude and search set radius')
disp([norm(K),B(1)])
K'
% Compute the pid gains
kp=K(1)-K(2)*T;ki=1/K(2);kd=K(3)-K(1)*T+K(2)*T^2;
if pos==0;ki=ki*kp;kd=kd/kp;end

npid=K(1)*npid1+K(2)*npid2+K(3)*npid3; 
disp('PID zeros')
disp(roots(npid))

while length(npid) < length(dpid),npid=[0 npid];end
PID=[npid;dpid;kp,ki,kd];
if fplantf==0
  while length(numo) < length(deno),numo=[0 numo];end
end

if flag==1
clg;
loglog(frw(1,:),abs(WK));
title('Approximation Error');pause
x=pid_eval(numo,deno,PID,FIL,lam,frw(1,:),fplantf);
end

if B(5)>0
lx=lam*1.1;
% FIL=fil_tune(lam,al,lx,numo,deno,PID,B(5),frw);

if flag==1
% x=pid_eval(numo,deno,PID,FIL,lam,frw);
end
end
