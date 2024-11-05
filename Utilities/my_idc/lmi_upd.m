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
