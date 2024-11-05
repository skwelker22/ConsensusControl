function [theta,P,R]=adils(theta,P,R,y,w,gam,alp,Q0,thresh,thl,thh,pdz);

%function [theta,P]=adils(theta,P,R,y,w,gam,alp,Q,thresh,thl,thh,pdz);
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
%  pdz: projection dead-zone; 
%       active when pdz>0, projects so that the a posteriori
%       error is below pdz; max iterations can be defined by pdz(2);
%uses: lsproj
%Outputs: 	theta,P

% initializations
nth=length(theta);ny=length(y);
It=eye(nth,nth);I1=eye(ny,ny);
max_iter=10;toler=1e-4;

if nargin < 5
gam=0;alp=0;Q0=0;thresh=0;thl=[];thh=[];pdz=0;
elseif nargin < 6
alp=0;Q0=0;thresh=0;thl=[];thh=[];pdz=0;
elseif nargin < 7
Q0=0;thresh=0;thl=[];thh=[];pdz=0;
elseif nargin < 8
thresh=0;thl=[];thh=[];pdz=0;
elseif nargin < 9
thl=[];thh=[];pdz=0;
elseif nargin < 11
pdz=0;
end

   if gam == 0 ; gam=1; end
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
   if isempty(pdz) | proji==0;pdz=0;end

if length(pdz)>1;max_iter=pdz(2);end

% parameter updates
P1=gam*P+2*w*w'; R1=gam*R+2*y*w;
S=R1-P1*theta;Pin=pinv(P1);
e1=S'*Pin*S/max(1,trace(P1)); indexp=0;
if abs(e1) > thresh
   P=P1;R=R1;
   theta=theta+Pin*S;
% projection
      if proji>0
         [theta,index]=lsproj2(theta,thhi,thlo,P,Pin,toler);
         if indexp<index;indexp=index;end
      end
end

