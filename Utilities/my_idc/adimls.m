function [theta,P,R,J,M,mmeigH]=adimls(theta,P,R,J,M,y,w,m,lam,qtol,dz,thl,thh,pdz);

%function [theta,P,R,J,M]=adimls(theta,P,R,J,M,y,w,m,lam,qtol,dz,thl,thh,pdz);
%Least-Squares with fading-memory adaptation (discrete-time)
%
%Inputs:	
%		theta, P, R: parameter estimates, covariance, cross cov 
%		J:	cost functional
%		y, w: 	output and	regressor
%Optional inputs: (enter zero to select default values)
%       m:  normalization signal input and floor
%		lam: 	adaptation gain (default = 1)
%		qtol: covariance inversion floor (default = pinv)
%		dz: adaptation error threshold (default = 0)
%		thl,thh: parameter bounds (low-high)
%          ls-projection is active when supplied
%  pdz: projection dead-zone; 
%uses: lsproj
%Outputs: 	theta,P

% initializations
nth=length(theta);ny=length(y);
It=eye(nth,nth);I1=eye(ny,ny);
max_iter=10;toler=1e-4;aopt=-1;

if nargin < 13;thh=[];end
if nargin < 12;thl=[];end
if nargin < 11;dz=[];end;	if isempty(dz);dz=0;end
if nargin < 10;qtol=[];end
if nargin < 9;lam=[];end;	if isempty(lam);lam=1;end
if nargin < 8;m=[];end		

if isempty(m);m_u=0;m_off=1;
elseif length(m)<ny+1;m_u=m(1:ny);m_off=1;
else;m_off=m(ny+1);m_u=m(1:ny);
end

if isempty(thl) | isempty(thh);
    proji=0;
elseif length(thl)==1 & length(thh)==1
    proji=1;thlo=0*theta+thl;thhi=0*theta+thh;
elseif length(thl)==1 & size(thh)==size(theta)
    proji=1;thlo=0*theta+thl;thhi=thh;
elseif size(thl)==size(theta) & length(thh)==1
    proji=1;thhi=0*theta+thh;thlo=thl;
elseif size(thl)==size(theta) & size(thh)==size(theta)
    proji=1;thlo=thl;thhi=thh;
else
    proji=0;
end

P1=P;R1=R;S=0*R;DZC=zeros(ny,1);
% parameter updates
for i=1:ny
    p_ind=[nth*(i-1)+1:nth*i];
    R1(:,i)=lam*R(:,i)+2*y(i)*w(i,:)';
    P1i=lam*P(:,p_ind)+2*w(i,:)'*w(i,:);
    P1(:,p_ind)=P1i;
    S(:,i)=R1(:,i)-P1i*theta;
    if isempty(qtol);IP1i=pinv(P1i);else;IP1i=pinv(P1i,qtol);end
    DZC(i)=S(:,i)'*IP1i*S(:,i)/2;
end
J=lam*J+(y-w*theta).^2;
M=lam*M+m_u.*m_u;
M_m=max(M,m_off);
J_M=J./M_m; max_j=max(J_M);
a_i=zeros(ny,1);
ind_j=find(max_j-J_M<=0.1*max_j);
a_i(ind_j)=1+a_i(ind_j);
del_J=-S*(a_i./M_m);
del2_J=0*theta*theta';
for i=1:ny
    p_ind=[nth*(i-1)+1:nth*i];
    del2_J=del2_J+a_i(i)*P(:,p_ind)/M_m(i);
end
e1=DZC./M_m; indexp=0;
if max(e1) >= dz*dz
    P=P1;R=R1;
    if isempty(qtol);
        del2_Ji=pinv(del2_J);
    else
        del2_Ji=pinv(del2_J,qtol);
    end
    dth=real(-del2_Ji*del_J);
    theta1=theta+dth;
    % projection
    if proji>0
        [theta1,index]=lsproj2(theta1,thhi,thlo,del2_J,del2_Ji,toler);
        if indexp<index;indexp=index;end
    end
    dth=real(theta1)-theta;
    JJ0=J;JJ1=(-S'*dth);JJ2=0*J;
    for i=1:ny
        p_ind=[nth*(i-1)+1:nth*i]; P1i=P(:,p_ind);
        JJ2(i)=1/2*dth'*P1i*dth;
    end
    
    if ny>1
        a=[0:.01:1];
        JX=(JJ0./(max(M,m_off)))*(0*a+1)+...
            (JJ1./(max(M,m_off)))*a+...
            (JJ2./(max(M,m_off)))*(a.*a);
        JXmax=max(JX);
        ind_a=find(JXmax==min(JXmax));ind_a=max(ind_a);
        aopt=a(ind_a);
        %  aopt=1;
    else
        aopt=1;
    end
    
    theta=theta+aopt*dth;
    J=JJ0+aopt*JJ1+aopt*aopt*JJ2;
end

eigH=log(sqrt(max(eig(del2_J),1e-9)))/log(10);
mmeigH=([min(eigH);max(eigH);aopt]);
