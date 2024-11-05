%  Script file TPIDF
%
%  Given square plant MIMO model, performs
%    diagonal approximation and PID design via loop shaping
%  Some Early Choices:
%  Default frequency range .0005-100 
%  Filter PID with 1/[Ts+1];T=1(min)

   T_max=input('Simulation time [10]   ');
   if length(T_max)<1;T_max=10;end;
   T_der=input('Derivative filter TC  [0.01]   ');
   if length(T_der)<1;T_der=.01;end;
   frwx=logspace(log10(5.e-4),log10(100),500);


innersys=input('Give Inner loop system file (e.g., FIDs...)  ','s');
sch_tol=input('Give reduction tolerance  [1.e-6]');
if length(sch_tol)<1;sch_tol=1.e-6;end

eval(['load ' innersys])
if T_max<0;T_max=abs(T_max);Bpr=-Bpr;Dpr=-Dpr;end
SI_=mksys(Apr,Bpr,Cpr,Dpr,'ss');
MUNC1=MUNC;

[Noutp,Ninp]=size(Dpr);
if Noutp~=Ninp
disp('System not square: resetting I/O')
NNX=min(Noutp,Ninp);Noutp=NNX;Ninp=NNX;
Bpr=Bpr(:,1:NNX);Cpr=Cpr(1:NNX,:);Dpr=Dpr(1:NNX,1:NNX);
end


disp('Step 1:  extract inner diagonals and perform order reduction')

if Noutp >1
     IX=1;
     [Aprx,Bprx,Cprx,Dprx,tot,hsv]=...
        schmr(Apr,Bpr(:,IX),Cpr(IX,:),Dpr(IX,IX),2,sch_tol);
     IX=2;
     [Aprx2,Bprx2,Cprx2,Dprx2,tot,hsv]=...
        schmr(Apr,Bpr(:,IX),Cpr(IX,:),Dpr(IX,IX),2,sch_tol);
   [Apn,Bpn,Cpn,Dpn]=append(Aprx,Bprx,Cprx,Dprx,...
                         Aprx2,Bprx2,Cprx2,Dprx2);
else
   Apn=Apr;Bpn=Bpr;Cpn=Cpr;Dpn=Dpr;
end
if Noutp > 2
   for IX=3:Noutp
     [Aprx,Bprx,Cprx,Dprx,tot,hsv]=...
        schmr(Apr,Bpr(:,IX),Cpr(IX,:),Dpr(IX,IX),2,sch_tol);
     [Apn,Bpn,Cpn,Dpn]=append(Apn,Bpn,Cpn,Dpn,...
                         Aprx,Bprx,Cprx,Dprx);
   end
end

[Ax,Bx,Cx,Dx]=addss(Apr,Bpr,Cpr,Dpr,Apn,Bpn,-Cpn,-Dpn);
svPx=sv3_5(Ax,Bx,Cx,Dx,-1,frun);
svPxs=sv3_5(Ax,Bx,Cx,Dx,1,frun);svPxs=svPxs(1,:);
svPx=min([svPx;svPxs]);
svPP=sv3_5(Apn,Bpn,Cpn,Dpn,1,frun);
MUNCx=(svPx)'./svPP(Noutp,:)';
loglog(frun,sat([MUNCx,MUNC1(:,1)],1.e4,1.e-4),...
       frun,sat(MUNC1(:,1)+MUNCx,1.e4,1.e-4),'--');grid;
title('Approximation Error SV'); pause
%---------------------------------------------------------

disp('Step 2:  TUNE  PIDs ')
tc_filt=zeros(Noutp,1);KPIDIN=zeros(3*Noutp,3);tdone=0;
while tdone==0
lwdef=input('Give default value for freq.weight [1]  ');
if length(lwdef)<1;lwdef=1;end

for IX = 1:Noutp
disp('Tuning loop:'),disp(IX)

%------  MODEL REDUCTION --------
     [Aprx,Bprx,Cprx,Dprx,tot,hsv]=...
        schmr(Apr,Bpr(:,IX),Cpr(IX,:),Dpr(IX,IX),2,sch_tol);
%------  TUNING -----------------
     disp('roots of numerator and denominator')
     [NUMx,DENx]=ss2tf(Aprx,Bprx,Cprx,Dprx,1);
     roots(NUMx),roots(DENx)
     T_fil=input('Error Filter TC [0.1]   ');
     if length(T_fil)<1;T_fil=.1;end;pidfil=[T_fil 1];
     tc_filt(IX)=T_fil;
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
     KPIDIN((IX-1)*3+1:IX*3,:)=PIDx(1:3,:);
end
%--------  EVALUATING DESIGN  ------------------
kpidin=KPIDIN([3:3:3*Noutp],:);
    [Acr,Bcr,Ccr,Dcr]=npidss(kpidin,T_der,tc_filt,0);

%--------  EVALUATING DESIGN  ------------------

[Algn,Blgn,Clgn,Dlgn]=series(Acr,Bcr,Ccr,Dcr,Apr,Bpr,Cpr,Dpr);
[Acly,Bcly,Ccly,Dcly]=feedbk(Algn,Blgn,Clgn,Dlgn,2);

clg;TIM=[0:T_max/500:T_max]';YIS=zeros(length(TIM),Noutp*Noutp);
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
end

kpidin=KPIDIN([3:3:3*Noutp],:);


savsys=input('Give filename to save results (e.g., out...)  ','s');
if length(savsys)>2
  eval(['save ',savsys,' kpidin tc_filt T_der'])
end




