function [theta,P]=adfmls(theta,PP,y,w,gam,alp,Q0,thresh,thl,thh,pdz);

%function [theta,P]=adfmls(theta,P,y,w,gam,alp,Q,thresh,thl,thh,pdz);
%Least-Squares with fading-memory adaptation (discrete-time)
%
%Inputs:	
%		theta, P: 	parameter estimates and covariance 
%		y, w: 	output and	regressor
%Optional inputs: (enter zero to select default values)
%		gam: 	adaptation gain (default = 1)
%		alp: 	forgetting factor (0=gradient/newton,1=LS, def=0)
%		Q: 	covariance floor (default = eye(P)*1.e-4)
%		thresh: adaptation error threshold (default = 0)
%		thl,thh: parameter bounds (low-high)
%          ls-projection is active when supplied
%  pdz[2]: projection dead-zone; 
%       active when pdz(1)>0, projects so that the a posteriori
%       error is below pdz(1); max iterations are defined by pdz(2);
%uses: lsproj
%Outputs: 	theta,P

% initializations
nth=length(theta);ny=length(y);
It=eye(nth,nth);I1=eye(ny,ny);
max_iter=10;toler=1e-4;

if nargin < 5
    gam=0;alp=0;Q0=0;thresh=0;thl=[];thh=[];pdz=[0,max_iter];
elseif nargin < 6
    alp=0;Q0=0;thresh=0;thl=[];thh=[];pdz=[0,max_iter];
elseif nargin < 7
    Q0=0;thresh=0;thl=[];thh=[];pdz=[0,max_iter];
elseif nargin < 8
    thresh=0;thl=[];thh=[];pdz=[0,max_iter];
elseif nargin < 9
    thl=[];thh=[];pdz=[0,max_iter];
elseif nargin < 11
    pdz=[0,max_iter];
end

if gam < 0 ; gam=1; end
if Q0 == 0 ; Q=0.0001*It;else;Q=Q0*It; end
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
if length(pdz)==0;
    PDZ=[0,max_iter];
elseif length(pdz)==1; 
    PDZ=[pdz(1),max_iter];
else
    PDZ=[pdz(1),pdz(2)];
end
max_iter=PDZ(2);
if isempty(pdz) | proji==0;PDZ(1)=0;end

if PP == 0; P=It;else;P=PP*It;end

% parameter updates
e1=y-w'*theta; indexp=0;
if abs(e1) > thresh
    PIN=P\It; L1=1/(I1+abs(gam*w'*PIN*w)); e1k=L1*e1;
    theta=theta+gam*PIN*w*e1k;thetau=theta;
    % projection along constant error manifolds
    if PDZ(1)>0
        e0=y-w'*theta;
        [theta,index]=lsproj2(theta,thhi,thlo,P,PIN,toler);
        e1=y-w'*theta; tempe=e0;niter=0;
        IND=abs(e1);
        if indexp<index;index=indexp;end 
        while IND>PDZ(1)+abs(e0) & niter<max_iter
            niter=niter+1;
            theta=theta+gam*PIN*w*L1*e1;
            [theta,index]=lsproj2(theta,thhi,thlo,P,PIN,toler);
            if indexp<index; indexp=index;end
            e1=y-w'*theta;IND=abs(e1);
        end
    end
    % projection
    if proji>0
        [theta,index]=lsproj2(theta,thhi,thlo,P,PIN,toler);
        if indexp<index;indexp=index;end
    end
    % covariance updates 
    if indexp == 0
        P=alp*P+(1-alp)*Q+gam*alp*w*w';
    elseif abs(w'*(theta-thetau)-e1k) <= abs(e1k)
        P=alp*P+(1-alp)*Q+gam*alp*w*w';
    end
end

