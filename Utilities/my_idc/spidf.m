%  Script file SPIDF
%
%  Given square plant MIMO model, performs
%    diagonal approximation and PID design via loop shaping
%  Some Early Choices:
%  Default frequency range .0005-100 
%  Filter PID with 1/[Ts+1];T=1(min)

   T_max=input('Simulation time [10]   ');
   if length(T_max)<1;T_max=10;end;
   T_der=input('Derivative filter TC [0.01]  ');
   if length(T_der)<1;T_der=.01;end;
   frwx=logspace(log10(5.e-4),log10(100),200);


innersys=input('Give Inner loop system file (e.g., FIDs...)  ','s');

eval(['load ' innersys])
SI_=mksys(Apr,Bpr,Cpr,Dpr,'ss');
MUNC1=MUNC;

[Noutp,Ninp]=size(Dpr);
if Noutp~=Ninp
disp('Cannot handle non-square systems, yet')
end

disp('Step 1:  extract inner diagonals and perform order reduction')

[Apn,Bpn,Cpn,Dpn]=append(Apr,Bpr(:,1),Cpr(1,:),Dpr(1,1),...
                         Apr,Bpr(:,2),Cpr(2,:),Dpr(2,2));
if Noutp > 2
   for IX=3:Noutp
     [Apn,Bpn,Cpn,Dpn]=append(Apn,Bpn,Cpn,Dpn,...
                         Apr,Bpr(:,IX),Cpr(IX,:),Dpr(IX,IX));
   end
end

[Ax,Bx,Cx,Dx]=addss(Apr,Bpr,Cpr,Dpr,Apn,Bpn,-Cpn,-Dpn);
svPx=sv3_5(Ax,Bx,Cx,Dx,-1,frun);
svPxs=sv3_5(Ax,Bx,Cx,Dx,1,frun);svPxs=max(svPxs);
svPx=min([svPx;svPxs]);
svPP=sv3_5(Apn,Bpn,Cpn,Dpn,1,frun);
MUNCx=(svPx)'./min(svPP)';
loglog(frun,sat([MUNCx,MUNC1(:,1)],1.e4),...
       frun,sat(MUNC1(:,1)+MUNCx,1.e4),'--');grid;
title('Approximation Error SV'); pause
%---------------------------------------------------------

disp('Step 2:  TUNE  PIDs ')
tc_filt=zeros(Noutp,1);KPIDIN=zeros(3*Noutp,3);tdone=0;
while tdone==0
L_seq=input('give loop closing sequence  ');
lwdef=input('Give default value for freq.weight [1]  ');
if length(lwdef)<1;lwdef=1;end
IX=1;
disp('Tuning loop:'),disp(L_seq(IX))
     JX=L_seq(IX);
%------  MODEL REDUCTION --------
     max(real(eig(Apr)))
     delred=input('give delta-shift for mod-ord red  ')
       if length(delred)<1;delred=0;end
     Aprd=Apr-delred*eye(length(Apr),length(Apr));
     [Aprx,Bprx,Cprx,Dprx,tot,hsv]=...
        schmr(Aprd,Bpr(:,JX),Cpr(JX,:),Dpr(JX,JX),3);
     Aprx=Aprx+delred*eye(length(Aprx),length(Aprx));
%------  TUNING -----------------
     disp('roots of numerator and denominator')
     [NUMx,DENx]=ss2tf(Aprx,Bprx,Cprx,Dprx,1);
     roots(NUMx),roots(DENx)
     T_fil=input('Error Filter TC [0.1]  ');
     if length(T_fil)<1;T_fil=.1;end;pidfil=[T_fil 1];
     tc_filt(JX)=T_fil;

     lam1=input('give bandwidth and type of PID design  ');
     op_t=input('give optimization method 0=Hinf, val>1=H2/Hinf (def=0)  ');
        if length(op_t)<1;op_t=0;end
        if length(lam1)<2;lam1(2)=0;end
        if length(lam1)<3;lam1(3)=0.3;end
        if length(lam1)<4;lam1(4)=0;end
        if length(lam1)<5;lam1(5)=lwdef;end
     PIDx=PIDftune(lam1(1),NUMx,DENx,lam1(2:4),op_t,...
                   [T_der,T_fil],frwx,[0,-1,-1],lwdef,1);
     disp('PID gains')
     PIDx
     KPIDIN((JX-1)*3+1:JX*3,:)=PIDx(1:3,:);

    [Acr,Bcr1,Ccr1,Dcr1]=tf2ss(KPIDIN((JX-1)*3+1,:),...
               conv(KPIDIN((JX-1)*3+2,:),pidfil));
    Bcr=Bcr1*zeros(1,Noutp);Ccr=zeros(Noutp,1)*Ccr1;Dcr=Ccr*Bcr;
    Bcr(:,JX)=Bcr1;Ccr(JX,:)=Ccr1;Dcr(JX,JX)=Dcr1;
    L_index=eye(Noutp,Noutp);L_index(JX,JX)=0;
%--------------------------------------------------------------
for IX = 2:Noutp
disp('Tuning loop:'),disp(L_seq(IX))
     JX=L_seq(IX);
     Q=inv(eye(Noutp,Noutp)+Dcr*Dpr);
     Atot=[Apr-Bpr*Q*Dcr*Cpr,          Bpr*Q*Ccr;...
           Bcr*Dpr*Q*Dcr*Cpr-Bcr*Cpr,  Acr-Bcr*Dpr*Q*Ccr];
     Btot=[Bpr*Q*L_index;Bcr*Dpr*Q*L_index];
     Ctot=[Cpr-Dpr*Q*Dcr*Cpr,Dpr*Q*Ccr];
     Dtot=[Dpr*Q*L_index];

%------  MODEL REDUCTION --------
     max(real(eig(Atot)))
     delred=input('give delta-shift for mod-ord red  ')
       if length(delred)<1;delred=0;end
     Atotd=Atot-delred*eye(length(Atot),length(Atot));
     [Aprx,Bprx,Cprx,Dprx,tot,hsv]=...
        schmr(Atotd,Btot(:,JX),Ctot(JX,:),Dtot(JX,JX),3);
     Aprx=Aprx+delred*eye(length(Aprx),length(Aprx));
%------  TUNING -----------------
     disp('roots of numerator and denominator')
     [NUMx,DENx]=ss2tf(Aprx,Bprx,Cprx,Dprx,1);
     roots(NUMx),roots(DENx)
     T_fil=input('Error Filter TC [0.1]  ');
     if length(T_fil)<1;T_fil=.1;end;pidfil=[T_fil 1]
     tc_filt(JX)=T_fil;
     lam1=input('give bandwidth and type of PID design  ');
     op_t=input('give optimization method 0=Hinf, val>1=H2/Hinf (def=0)  ');
        if length(op_t)<1;op_t=0;end
        if length(lam1)<2;lam1(2)=0;end
        if length(lam1)<3;lam1(3)=0.3;end
        if length(lam1)<4;lam1(4)=0;end
        if length(lam1)<5;lam1(5)=lwdef;end
     PIDx=PIDftune(lam1(1),NUMx,DENx,lam1(2:4),op_t,...
                   [T_der,T_fil],frwx,[0,-1,-1],lwdef,1);
     disp('PID gains')
     PIDx
     KPIDIN((JX-1)*3+1:JX*3,:)=PIDx(1:3,:);
    [APID,BPID1,CPID1,DPID1]=tf2ss(KPIDIN((JX-1)*3+1,:),...
               conv(KPIDIN((JX-1)*3+2,:),pidfil));
    BPID=BPID1*zeros(1,Noutp);CPID=zeros(Noutp,1)*CPID1;DPID=CPID*BPID;
    BPID(:,JX)=BPID1;CPID(JX,:)=CPID1;DPID(JX,JX)=DPID1;
    [Acr,Bcr,Ccr,Dcr]=addss(Acr,Bcr,Ccr,Dcr,APID,BPID,CPID,DPID);
    L_index(JX,JX)=0;
end
%--------  EVALUATING DESIGN  ------------------

[Algn,Blgn,Clgn,Dlgn]=series(Acr,Bcr,Ccr,Dcr,Apr,Bpr,Cpr,Dpr);
[Acly,Bcly,Ccly,Dcly]=feedbk(Algn,Blgn,Clgn,Dlgn,2);

clf;TIM=[0:T_max/500:T_max]';YIS=zeros(length(TIM),Noutp*Noutp);
   for IX=1:Noutp
     YIS(:,(IX-1)*Noutp+1:IX*Noutp)=step(Acly,Bcly,Ccly,Dcly,IX,TIM);
   end
plot(TIM,YIS);title('Step Responses');pause
Y=lsim(Acly,Bcly,Ccly,Dcly,ones(length(TIM),Noutp),TIM);plot(TIM,Y);
title('All-channel Step Response');pause
[LPMG] = abs(freqrc(Acly,Bcly,Ccly,Dcly,frun));
loglog(frun,LPMG);title('Frequency Responses');pause
%------------------------------------------------
tdone=input('Done ?   ');
if isempty(tdone);tdone=1;end
end

kpidin=KPIDIN([3:3:3*Noutp],:);


savsys=input('Give filename to save results (e.g., out...)  ','s');
if length(savsys)>2
  eval(['save ',savsys,' kpidin tc_filt T_der'])
end



