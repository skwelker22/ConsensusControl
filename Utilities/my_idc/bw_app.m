function cl_spec=bw_app(frun,MUNC,S_roll,RISK,AGRS,T_rollm,dis_f);
% function cl_spec=bw_app(frun,MUNC,S_roll,RISK,AGRS,T_rollm,dis_f);
% provides first-cut rule for the selection of the closed-loop
% specs based on the multiplicative uncertainty in MUNC

if nargin<3;S_roll=2;end
if isempty(S_roll);S_roll=2;end
if nargin<4;RISK=1;end
if nargin<5;AGRS=1;end
if nargin<6;T_rollm=3;end
if nargin<7;dis_f=1;end

if RISK<0; BWR=-RISK;RISK=1;else;BWR=0;end

n1=max(find(MUNC(:,1)<1));if isempty(n1);n1=1;end
n0=max(find(MUNC(:,1)<.2));
  if isempty(n0);n0=max(find(MUNC(:,1)<.5));end
  if isempty(n0);n0=max(find(frun<=0.1*frun(n1)));end
  if isempty(n0);n0=1;end
n2=min(find(MUNC(:,1)>10));
  if isempty(n2);n2=max(find(frun<=10*frun(n1)));end
  if n2-n1<5;n2=max(find(frun<=10*frun(n1)));end
n3=min(find(MUNC(:,1)>100));
  if isempty(n3);n3=min(find(MUNC(:,1)>40));end
  if isempty(n3);n3=n2;end
  if n3-n1<10;n3=max(find(frun<=30*frun(n1)));end
n4=length(frun);
n0x=max(find(frun<=0.1*BWR));if isempty(n0x);n0x=1;end
if BWR>0;n0=min(n0,n0x);end
n0y=max(find(frun<=10*BWR));
if BWR>0;n3=max(n3,n0y);end

lfr=log(frun(n1:n4));lmunc=log(MUNC(n1:n4,1));
A=[lfr,0*lfr+1];x=A\lmunc;
slope1=round(abs(x(1))+.3);

lfr=log(frun(n0:n2));lmunc=log(MUNC(n0:n2,1));
A=[lfr,0*lfr+1];x=A\lmunc;
slope2=round(abs(x(1))+.3);
%disp([slope1,slope2])
slope=max([slope1,slope2,1]);
slope=abs(min(slope,T_rollm));

A=[0*lfr+1];B=lmunc-slope*lfr;x=A\B;
BW_nom=exp(-x/slope);BW=BW_nom*.7;
if BWR>0; BW_nom=BWR;end

ssTh=awgt_sel(frun,BW_nom,slope,0);
ssT=awgt_sel(frun,BW,slope,0);
err=(max(ssT(n0:n3).*MUNC(n0:n3,1))-RISK*0.3);
kiter=0;
while abs(err)>.01
   kiter=kiter+1;
   BW=exp(log(BW)-err/(1+abs(err)));
   if BWR>0;BW=BWR;end
   ssT=awgt_sel(frun,BW,slope,0);
   err=(max(ssT(n0:n3).*MUNC(n0:n3,1))-RISK*0.3);
   if kiter > 60;err=0;disp('bw_app: max. iterations exceeded');end
   if BWR>0;err=0;end
end


cl_spec=[BW,slope,0,S_roll,RISK,AGRS];

sbw=sbw_rule(cl_spec);
cl_spec(3)=sbw;

