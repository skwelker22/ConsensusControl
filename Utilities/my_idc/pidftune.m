function [PID,Ub]=PIDFtune(lam,num,den,lp_s,op_t,f_,fr_r,con_,fr_w,flag);
%  Usage:
%  [PID,Ub]=PIDFtune(lam,numo,deno,[lp_s,op_t,f_,fr_r,con_,fr_w,flag]);
%
%  Computes the PID gains kp,ki,kd, and a filter 1/(TF*s+1), 
%      tuned  by loop-shaping techniques to achieve a closed-loop
%      bandwidth approx. equal to lam,
%      while maintaining a 'small' sensitivity peak value.
%
%  Limitations: The plant must be stabilizable by a PID with 
%      positive gains. Also, suitable choice of the target loop 
%      may be required for plants with integral action to achieve
%      adequate input disturbance attenuation.
%
% PID definition: 
%  Factorized form:  u = [kp*(1+1/[ki*s]+kd*s/[T*s+1])]*[1/(TF*s+1)] e
%  Additive form: u = [kp+1/[ki*s]+kd*s/[T*s+1]]*[1/(TF*s+1)] e
%
% Inputs:  
%   lam: desired bandwidth; constrained by robustness
%        spec's and modelling errors; 
%   num,den: numerator,denominator of nominal plant.
%        The plant can also be specified by its frequency response;
%        in this case, set num=0 and den=row vector of freq.resp.
%        evaluated at the frequency points in fr_r (e.g., by freqs);
%        also, for this case, the option lp_s(1)=2 cannot be used. 
%   lp_s:  Row vector specifying the target loop shape:
%      lp_s(1)=loop shape type, 
%         =0 selects a single integrator as a target loop 
%            (OK for plants with no poles near the jw-axis)
%         ~=0 selects higher order target loop with pole/zero 
%            specified by  according to the following convention:
%                1-3: L(s)=lam(s+a)/s(s+p), 
%                  4: L(jw) in fr_r(2,:); (row vector)
%            The values of a and p are selected according to the table 
%         =1: a=lam*lp_s(2), p=0; 
%         =2: a=lam*lp_s(2), p are the smallest in magnitude plant pole. 
%           (1+2 appropriate for plants with a pole near the jw-axis;
%            here larger values of lp_s(2) may induce sensitivity 
%            peaking; smaller values may affect input disturbance 
%            attenuation; typical choices are around 0.2-0.3)
%         =3: a=lp_s(2);p=lp_s(3) (general form)
%         =4: L(jw) in fr_r(2,:); in this case lam is the expected BW.
%         =5: L(jw) from lqr: k=lqr(A,B,Q,R,N); L=[A,B,K,0]
%             Q=C'C+lp_s(2)I, R=lp_s(3)I+D'D; if lp_s(3)=0, iterate until S-BW=lam
%   op_t: Row Vector of switches specifying the optimization type;
%         op_t(1)=0 performs H-inf optimization,
%                ~=0 (=val) performs H-2 opt. with the H-inf constraint
%                  ||(Sens.)(loop_t.f.-target)||_inf < val*min||.||_inf
%         op_t(2)=0/1 selects PID/PI controllers
%   f_:  Row Vector of filter parameters
%         f_(1)= filter time constant for the pseudo-derivative term;
%                typically 0.5 or 1 * sampling time; default = 0.01;
%         f_(2)= TF, the filter t.c. (default=0)
%   fr_r: Frequency range for the optimization;
%         =0: 300 pts in 0.01*lam--100*lam;
%         =1: 200 pts in 0.1*lam--10*lam;
%         =2: 500 pts in 0.001*lam--1000*lam;
%         or, it may contain a row vector of user specified frequencies;
%         if desired, its second row can be used to specify the 
%           freq. response of the loop shape at the same frequencies;
%         if the plant is given as a frequency response (den), 
%           fr_r should be the corresponding frequencies.
%   con_: Row vector specifying the constraints:
%         con_(1)= search set radius; a warning message is displayed
%                  if the optimizer exceeds this value;
%            NOTE: Larger values may slow down convergence, 
%                  smaller values may overly restrict the optimizer.
%         con_(2)= High frequency gain constraint (disabled by -1)
%         con_(3)= Derivative term gain constraint (disabled by -1)
%         con_(4)= 0/1 selects summation/factor form of the PID
%                  1 restricts the PID gains to be positive; 
%                  0 restricts the PID to be minimum phase only;
%         Default:[1.e3,1.e6,1.e6,0].
%             NOTE: H.F. gain (summation form)=(kp*T+kd)/(T*TF) or
%               (kp*T+kd)/T; Deriv. term gain (summation form)=kd;
%   fr_w: Frequency dependent fitting weight 
%        if scalar <1, low frequencies are given higher weight; 
%        if scalar >1, high frequencies are given higher weight;
%        a row vector [f1 f2] selects the filter fr_w(s)=(s+f1)/(s+f2); 
%        mid-frequencies can be emphasized with a row vector [f1,f2,f3]
%          fr_w(s)=(s+f1/f2)(s+f1*f2)/(s+f1/f3)(s+f1*f3); 
%          here, a first choice is [2lam, 2, 1]; NOTE:f2>f3;
%        the fitting weight can also be specified by using
%          fr_w(1,:)=num_weight;fr_w(2,:)=den_weight, or as a row
%          vector freq. response at the frequencies in fr_r; default=1;
%
%  Zero entries select default values; also incomplete vectors are
%  padded with zeros;
%
%  Intermediate Outputs:
%  - Ub: fitting error; must be less than 1 to guarantee closed loop 
%          stability; typical value should be around 0.2-0.3,
%          otherwise large sensitivity peaks may occur; 
%        if too high, a suitable scaling can be achieved by
%        reducing lam and kp (kp~lam~Ub); tuning should be 
%        repeated when lam is scaled by 1/2 order of magnitude.)
%  - PID zeros and Closed-loop poles
%  
%  Displays:
%       Closed loop sensitivity, complementary sensitivity,
%       control sensitivity and response to step reference
%       commands and disturbances
%
% OUTPUT PID format: [ numerator polynomial;
%                      denominator polynomial;
%                      kp, ki, kd;
%                      filter denominator
%                      out_flag, b, c (set-pt. weight)  ]
%
%    out_flag is set if the optimizer falls outside the search set radius


% K. TSAKALIS, 8/23/96

FIL=[1;1];al=[];pm=[];out_flag=0;
format short e

if nargin<4, lp_s=0; end
if nargin<5, op_t=0;end
if nargin<6, f_=0;end
if nargin<7, fr_r=0;end
if nargin<8, con_=0;end
if nargin<9, fr_w=0;end
if nargin<10, flag=1;end

if lp_s==0;lp_s=[0,0,0];end
if op_t==0;op_t=[0,0,0];end
if f_==0;f_=[0,0];end
if con_==0;con_=[0,0,0,0];end
if fr_w==0;fr_w=1;end
while length(op_t)<3;op_t=[op_t,0];end
while length(f_)<2;f_=[f_,0];end
while length(con_)<4;con_=[con_,0];end
while length(lp_s)<3;lp_s=[lp_s,0];end
f_=abs(f_);con_(1)=abs(con_(1));op_t(1)=abs(op_t(1));

if lp_s(1) == 1 & lp_s(2) <=0; lp_s(2)=.3;end
if lp_s(1) == 2 & num == 0 ; lp_s(1)=0;end
[nfrr,mfrr]=size(fr_r);
  if lp_s(1) == 4 & nfrr ~= 2 lp_s(1)=0;
     disp('PIDFTUNE WARNING: Loop shape specification error;') 
     disp('                  resetting to default')
  end 
if con_(1) == 0;con_(1)=1.e3;end
if con_(2) == 0;con_(2)=1.e6;end
if con_(3) == 0;con_(3)=1.e6;end

if f_(1)==0;f_(1)=0.001;end;T=f_(1);
TF=f_(2);

if fr_r == 0
   fmin=0.01*lam;fmax=100*lam;
   frw=logspace(log10(fmin),log10(fmax),300);
elseif fr_r==1
   fmin=0.1*lam;fmax=10*lam;
   frw=logspace(log10(fmin),log10(fmax),200);
elseif fr_r==2
   fmin=0.001*lam;fmax=1000*lam;
   frw=logspace(log10(fmin),log10(fmax),500);
end
if length(fr_r)>2;frw=fr_r(1,:);end

[nnls1,nnls2]=size(fr_r);
if nnls1 >1 lp_s(1)=4;end
if num==0; fplantf=1; else, fplantf=0;end
if fplantf ==0
  r2=roots(den);pmin=r2(min(find(abs(r2)==min(abs(r2)))));
  if real(pmin)==0;pmin=abs(pmin);else;pmin=sign(real(pmin))*abs(pmin);end
  r2=roots(num);
  if length(r2)>0;
     zmin=r2(min(find(abs(r2)==min(abs(r2)))));zmin=abs(zmin(1));
  else, 
     zmin=[];
  end
end

if lp_s(1)==1,al=lp_s(2)*lam;pm=0;end
if lp_s(1)==2,pm=-pmin;al=lp_s(2)*abs(lam+2*(pm));end
if lp_s(1)==3,al=lp_s(2);pm=lp_s(3);end
if lp_s(1) == 0
   lshape=freqs(lam,[1 0],frw).'; 
elseif lp_s(1)~= 4
        lshape=freqs(lam*[1/al 1],[1 pm 0],frw).';
else
   lshape=(fr_r(2,:)).';
end
if lp_s(1) == 5
   [n_L,d_L]=lqr_lsh(num,den,lam,lp_s);
   lshape=freqs(n_L,d_L,frw).'; 
end
% create frequency responses
     if fplantf==0
       pla=freqs(num,den,frw).';
     else
       pla=den.';
     end

if flag >= 1
disp('Loop shape parameters lam -zero -pole')
disp([lam,al,pm])
end

[nlw,mlw]=size(fr_w);
     if length(fr_w)==1
       fr_w=abs(fr_w);
       if fr_w<1
        lwwei=(freqs([1 lam],[1 fr_w*lam],frw)).';
       else
        lwwei=(freqs([1/lam 1],[1/(lam*fr_w) 1],frw)).';
       end
     elseif nlw==2
       lwwei=(freqs(fr_w(1,:),fr_w(2,:),frw)).';
     elseif mlw==2
        lwwei=(freqs([1 fr_w(1)],[1 fr_w(2)],frw)).';
        lwwei=lwwei/min(abs(lwwei));
     elseif mlw==3
       lwnum=conv([1,fr_w(1)/fr_w(2)],[1,fr_w(1)*frw(2)]);
       lwden=conv([1,fr_w(1)/fr_w(3)],[1,fr_w(1)*frw(3)]);
       lwwei=freqs(lwnum,lwden,frw).';
     else
       lwwei=fr_w.';
     end

%----- DEFINE FILTER ITERATION
  dpid=[T 1 0];npid1=[0 1 0];npid2=[0 0 1];npid3=[1 0 0];
  zp1=freqs(npid1,dpid,frw).';zp2=freqs(npid2,dpid,frw).';
  zp3=freqs(npid3,dpid,frw).'; 
  wei=lwwei./(1+lshape);zi=lshape.*wei;
  w1=zp1.*pla.*wei;w2=zp2.*pla.*wei;w3=zp3.*pla.*wei; 
  f_(2)=TF;  pidflt=freqs(1,[TF,1],frw).';
  W=[w1.*pidflt,w2.*pidflt,w3.*pidflt];
%-----------------------------------------------------
  op_ty=op_t;op_ty(1)=0;
  [K,Ub]=lmi_pidf(W,zi,2.e-5,f_,op_ty,con_,frw',flag);
  wei0=(1)./(1+lshape);zi0=wei0.*lshape; 
  W0=[zp1.*pla.*wei0.*pidflt,zp2.*pla.*wei0.*pidflt,...
      zp3.*pla.*wei0.*pidflt];
  WK=W0*K-zi0;

    if norm(K)>=con_(1)
     out_flag=1;
     if flag>=0
       disp('PIDFTUNE WARNING:')
       disp('Optimizer outside initial search set')
       disp('If desired, repeat with a higher value in con_(1)') 
       disp('Optimizer magnitude and search set radius')
       disp([norm(K),con_(1)])
     end
    end

% Compute the pid gains
kp=K(1)-K(2)*T;ki=K(2);kd=K(3)-K(1)*T+K(2)*T^2;
if con_(4)==1;ki=kp/ki;kd=kd/kp;end

npid=K(1)*npid1+K(2)*npid2+K(3)*npid3; 
   if flag >=1
     disp('PID zeros')
     disp(roots(npid))
   end

while length(npid) < length(dpid),npid=[0 npid];end
PID=[npid;dpid;kp,ki,kd;0,TF,1];
if fplantf==0
  while length(num) < length(den),num=[0 num];end
end

   if flag >=0
      figure(1);clf;loglog(frw(1,:),abs(WK));
      title('Approximation Error');pause
      disp('Approximation error and filter TC')
      disp([Ub,TF])
   end
   if flag >=1
      x=pidfeval(num,den,PID,TF,FIL,lam,frw,fplantf,flag);
   end

%--------------------------------------------------------
if op_t(1)~=0
  op_t(1)=max(op_t(1),1.001);
  op_ty(1)=-op_t(1)*Ub;
  [K,Ub]=lmi_pidf(W,zi,2.e-5,f_,op_ty,con_,frw',flag);
  wei0=(1)./(1+lshape);zi0=wei0.*lshape; 
  W0=[zp1.*pla.*wei0.*pidflt,zp2.*pla.*wei0.*pidflt,...
      zp3.*pla.*wei0.*pidflt];
  WK=W0*K-zi0;

    if norm(K)>=con_(1)
     out_flag=1;
     if flag>=0
       disp('PIDFTUNE WARNING:')
       disp('Optimizer outside initial search set')
       disp('If desired, repeat with a higher value in con_(1)') 
       disp('Optimizer magnitude and search set radius')
       disp([norm(K),con_(1)])
     end
    end

% Compute the pid gains
kp=K(1)-K(2)*T;ki=1/K(2);kd=K(3)-K(1)*T+K(2)*T^2;
if con_(4)==1;ki=kp/ki;kd=kd/kp;end

npid=K(1)*npid1+K(2)*npid2+K(3)*npid3; 
   if flag >=1
     disp('PID zeros')
     disp(roots(npid))
   end

while length(npid) < length(dpid),npid=[0 npid];end
PID=[npid;dpid;kp,ki,kd;0,TF,1];
if fplantf==0
  while length(num) < length(den),num=[0 num];end
end

   if flag >=0
      figure(1);clf;loglog(frw(1,:),abs(WK));
      title('Approximation Error');pause
      disp('Approximation error and filter TC')
      disp([Ub,TF])
   end
   if flag >=1
      x=pidfeval(num,den,PID,TF,FIL,lam,frw,fplantf,flag);
   end
end
%--------------------------------------------------


BC_=[0 1 1];

PID=[PID;BC_];
PID(5,1)=out_flag;
