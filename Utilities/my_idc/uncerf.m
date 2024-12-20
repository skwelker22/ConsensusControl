function [frun,AUN,MUNC,sc,rbw]=uncerf(A,F,B,C,D,By,Dy,X0,u,y,t,fpts,flag,cl_spec,ssW,ssDmi);

% [frun,AUNC,MUNC,sc,rbw]=uncerf(A,F,B,C,D,By,Dy,X0,u,y,t,fpts,flag,cl_spec,ssW,ssDmi);
% function UNCERF to compute an estimate of the 
% input multiplicative, and divisive uncertainty
% in sys.id. (MIMOID). Uses a simplified fft approach.
%    For general use, it requires the system in (nonminimal form)
%      [Apr,Bpr,Cpr,Dpr,Fpr,Bypr,Dypr];
%      i/o pair in [u,y], time index in t; X0 = IC
%      fpts = no. of points where the unc. is evaluated in fpoints.
%      cl_spec = closed loop spec parameters
%      ssW = output weight sv
%      ssDmi = inner loop multiplicative uncertainty (input weight) 
%    Outputs: [frun, AUNC,MUNC,sc,rbw]; MUNC=[mult., divis.] unc.
%      sc = estimated stability condition; rbw = resulting cl.lp. specs

if nargin<16;ssDmi=[];end
if nargin<15;ssTi=[];end
if nargin<14;cl_spec=[];end
if nargin<13;flag=1;end

stp=t(2)-t(1);rbw=cl_spec;
[nn,ninp]=size(u);[nn,noutp]=size(y);
nix=min(ninp,noutp);
i_nnan=find(~isnan(u(:,1)));
nan_ind=[find(isnan(u(:,1)));nn+1];
%-----------------------------------------------
% construct errors
%-----------------------------------------------
wu=lsim_nan(F,B-By*D,C,D,u,t,X0);
wy=lsim_nan(F,By,C,Dy,y,t);
Ferr=y-(wu+wy);
figure(2);clf,plot(t,Ferr);title('Estimation Error')
%-----------------------------------------------
disp('Uncertainty bound estimation ');
%-----------------------------------------------
% basic additive and multiplicative
%-----------------------------------------------

[frun,s_e,s_u,s_y]=unc_fft(Ferr,u,y,fpts,stp);
AUNC=s_e./s_u;       %unc_bn(s_e,s_u,s_u,1,1);
AUNCuraw=s_e./s_u;   %unc_bn(s_e,s_u,s_u,1,1);
AUNCyraw=s_e./s_y;   %unc_bn(s_e,s_y,s_y,1,1);
sspla=sv3_5(A,B,C,D,1,frun);
ssdf=sv3_5(F,By,-C,-Dy,3,frun);
ssnf=sv3_5(F,B-By*D,C,D,1,frun);
%AUNC=AUNC./(ssdf(noutp,:)');
MUNCo=AUNC./(sspla(nix,:)')./(ssdf(noutp,:)');
MUNCoLB = AUNC./(sspla(1,:)')./(ssdf(noutp,:)');
MUNCi=unc_bn(s_e,s_u,s_u,(1)./(ssnf(nix,:)),1);
FUNC=unc_bn(s_e,s_y,s_y,(1)./(ssdf(noutp,:)),1);

figure(5);clf
loglog(frun,MUNCi,'g',frun,MUNCoLB,'m',frun,MUNCo,'b',frun,FUNC,'r');grid;
title('Output (b)/LowerBound (m)/Input (g) Multiplicative and Divisive (r) Uncertainty')
%-----------------------------------------------
% Inner loop uncertainty
%-----------------------------------------------
if length(ssDmi)==length(frun)
    [nti,temp_n]=size(ssTi);nixt=min(nix,nti);
    %  inn_u=sspla(1,:)./sspla(nix,:).*(ssDmi');              % approximate
    %  inn_u=ssnf(1,:)./sspla(nix,:)./ssdf(nix,:).*(ssDmi');  % conservative 
    inn_u=ssnf(1,:)./ssnf(nix,:).*(ssDmi');                % T = scalar * I
    inn_f=(1+ssDmi');
else
    inn_u=0*frun';inn_f=inn_u+1;
end
MUNCoa=MUNCo.*(inn_f')+inn_u';

%figure(6);clf
%loglog(frun,MUNCo,'m',frun,sat(inn_u,1.e4,1.e-4),'y',frun,inn_f,'r',frun,MUNCoa,'b');
%grid;
%title('Output Multiplicative (m), Inner loop contribution (y,r) Total (b)')
%-----------------------------------------------
% Bandwidth estimation
%-----------------------------------------------
if length(cl_spec)>=4;S_roll=cl_spec(4);else;S_roll=2;end
if length(cl_spec)>=5;RISK=cl_spec(5);else;RISK=1;end
if length(cl_spec)>=6;AGRS=cl_spec(6);else;AGRS=.9;end
cl_speco=cl_spec;

if flag >0 | length(cl_spec)< 6
    cl_spec=q_spec(frun,MUNCoa,S_roll,RISK,AGRS);
else
    cl_specn=bw_app(frun,MUNCoa,S_roll,RISK,AGRS);
    rat=cl_specn(1)/cl_spec(1);
    if rat >1.25 | rat < .75
        disp('**  Large Deviation from Preset Specifications')
        disp(cl_speco),  disp(cl_specn) 
        if flag<0
            cl_spec=cl_speco;
        else
            cl_spec=q_spec(frun,MUNCoa,S_roll,RISK,AGRS);
        end
    end
end
%-----------------------------------------------
% almost necessary conditions
%-----------------------------------------------
done=1;
while done==1
    if flag>=1
        disp('Uncerf: current c-l specs [T_bw,T_roll,S_bw,S_roll,R,A]')
        disp(cl_spec)
        bwdth=input('specify target loop T,S ? (1=yes) [0]  ');
        if isempty(bwdth);bwdth=0;end
        if bwdth==1
            [j1,j2,j3,ft1,Tsens]=wgt_sel(1,frun,MUNCoa,0,flag);
            cl_spec(1)=min(j1);cl_spec(2)=length(j1);
            [j1,j2,j3,ft1,Sens]=wgt_sel(1,frun,[1.0./Tsens,MUNCoa,FUNC],-1,flag);
            cl_spec(3)=max(j1);cl_spec(4)=length(j1);
        else
            Tsens=awgt_sel(frun,cl_spec(1),cl_spec(2),0);
            Sens=awgt_sel(frun,cl_spec(3),cl_spec(4),-1);
        end
    else
        Tsens=awgt_sel(frun,cl_spec(1),cl_spec(2),0);
        Sens=awgt_sel(frun,cl_spec(3),cl_spec(4),-1);
    end
    Tsens=Tsens';Sens=Sens';
    %-----------------------------------------------------------
    
    ss1=Tsens./(ssnf(nix,:)).*inn_f;ss2=Sens./ssdf(noutp,:);
    ss11=Tsens./ssdf(noutp,:)./sspla(nix,:).*inn_f;
    
    figure(6);clf
    loglog(frun,Tsens,frun,Sens,frun,[1.0./MUNCoa 1./FUNC]);grid
    title('Weights and preliminary uncertainties')

    
    ssva=[ss1;ss2];
    UNCFi=unc_bn(s_e,s_u,s_y,ssva,2);
    MUNCi2=UNCFi(:,1)./(ssnf(nix,:)').*inn_f'+inn_u';
    FUNCi2=UNCFi(:,2)./(ssdf(noutp,:)');
    
    ssva=[ss11;ss2];
    UNCFo=unc_bn(s_e,s_u,s_y,ssva,2);
    MUNCo2=UNCFo(:,1)./(ssdf(noutp,:)'.*sspla(nix,:)').*inn_f'+inn_u';
    FUNCo2=UNCFo(:,2)./(ssdf(noutp,:)');
    
    %    loglog(frun,sat(MUNCo,1e5,1e-5),'b',...
    %        frun,sat(MUNCo2,1e5,1e-5),'r',...
    %        frun,sat(MUNCi,1e5,1e-5),'g',...
    %        frun,sat(MUNCi2,1e5,1e-5),'m');
    %    grid
    %    title('B,G=basic;R,M=weighted;B,R=output;G,M=input')
    %    if flag <=2;pause(2);else;pause;end
    %    loglog(frun,sat(FUNC,1e5,1e-5),'b',...
    %        frun,sat(FUNCo2,1e5,1e-5),'r',...
    %        frun,sat(FUNCi2,1e5,1e-5),'m');
    %    grid
    %    title('B,G=basic;R,M=weighted;B,R=output;G,M=input')
    %    if flag <=2;pause(2);else;pause;end
    
    %----------------------------------------------------
    R_A=inv_spec(cl_spec,frun,MUNCo2);
    if length(cl_spec)<6
        cl_spec(5:6)=R_A;
    else
        if abs(log(R_A(1)/cl_spec(5)))+abs(R_A(2)-cl_spec(6))>0.002
            cl_spec(5:6)=R_A;
        end
    end
    disp('Adjusted cl_spec')
    disp(cl_spec)
    %-----------------------------------------------------
    
    stabchk1=(Tsens'.*MUNCi2)+(Sens'.*FUNCi2); % --- needs upd
    stabchk2=(Tsens'.*MUNCo2)+(Sens'.*FUNCo2);
    if length(ssW)==length(frun)
        stabchk3=stabchk2*(max(ssW)./min(ssW))';
    else
        stabchk3=stabchk2;
    end
    
    figure(7);clf
    loglog(frun,stabchk1,frun,stabchk2,frun,stabchk3);grid
    title('stability checks')
    
    if flag>=1
        done=input('Repeat with different weights? (1=yes) [0]  ');
        if isempty(done);done=0;end
        if done==1;flag=2;end
    else
        done=0;
    end
end

if flag>=2
    isav=input('Saving estimates: 0,1=weighted o,i; 2,3 basic o,i  [0]  ');
    if length(isav)<1;isav=0;end
else
    isav=0;
end

if isav==0;MUNC=[MUNCo2,FUNCo2];AUN=[AUNC,UNCFo,AUNCuraw,AUNCyraw];end 
if isav==1;MUNC=[MUNCi2,FUNCi2];AUN=[AUNC,UNCFi,AUNCuraw,AUNCyraw];
    disp('Using Input Multiplicative Unc. requires sensitivity matching')
end 
if isav==2;MUNC=[MUNCo,FUNC];AUN=[AUNC,UNCFo,AUNCuraw,AUNCyraw];end 
if isav==3;MUNC=[MUNCi,FUNC];AUN=[AUNC,UNCFi,AUNCuraw,AUNCyraw];
    disp('Using Input Multiplicative Unc. requires sensitivity matching')
end 
sc=stabchk2;
rbw=cl_spec;
