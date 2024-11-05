function [K,Ub]=lmiqpid(W,zi,toler,con_,T,c_,flag);
% function [K,Ub]=LMIqPID(W,zi,toler,con,T,c_,flag);
%  Convex optimization function to perform PID tuning

% K. TSAKALIS, 8/23/01

flag =0;

K=[0;0;0]; Aell=con_(1)*con_(1)*eye(length(K),length(K));
L=[];U=[];err=1;erro=1;kount=0;phi=0;kkount=0;indp=0;

while err>toler
    [h,psi]=pidqcns(K,1.e-5,T,c_);
    h0=norm(h);
    if h0>0
        % 'constraint-iteration'
        iterty=-1;deepcut=psi;
    else
        % 'objective-iteration'
        iterty=1;
        [h,phi]=pidqobj(W,zi,K);
        erro=sqrt(h'*Aell*h);deepcut=0;
        %     L=[L;phi-erro];U=[U;phi];
        if length(U)<1;U=phi;else;U=min(U,phi);end
    end
    [K,Aell,err_flag]=lmi_upd(K,h,Aell,deepcut,flag);
    kount=kount+1;kkount=kkount+1;
    h0=norm(pidqcns(K,1.e-5,T,c_));
    err=(h0)+erro;
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
if K(1)<=tol,h(1)=-1;val=tol-K(1);
elseif K(2)<=tol,h(2)=-1;val=tol-K(2);
elseif K(3)<=tol,h(3)=-1;val=tol-K(3);
end
if K(3)-T*K(1)+T*T*K(2) >= cc
    h=[-T, T*T,1]';
    val=(K(3)-T*K(1)+T*T*K(2))-cc;
end

function [h,phi]=pidqobj(W,zi,K);
% function [h,phi]=PIDqOBJ(W,zi,K);
% Defines the optimization objective in pidftune

% K. TSAKALIS, 8/11/04


temp = W*K-zi;npt=length(temp);

phi  = max(abs(temp))^2;
indp = min(find(abs(temp)==sqrt(phi)));
h    = 2*real(temp(indp)'*W(indp,:));
h    = h'/norm(h);

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
