%
%  Given square plant MIMO models, performs
%    diagonal approximation and PID design via loop shaping
%  Default frequency range 5e-4-300 
%  Pseudo-derivative filter: s/[T_der s +1]; T_der=0.1
%  Filter PID with 1/[Ts+1];T=0.1(min)
%  PID gains are saved together with filter parametes.
%  Use "npidss" to obtain a state-space description

flag=1;
i_mult1=input('Input Multipliers [1]   ');
if isempty(i_mult1);i_mult1=1;end
T_max=input('Simulation time [10]   ');
if length(T_max)<1;T_max=10;end;
T_der=input('Derivative filter TC  [0.1]   ');
if length(T_der)<1;T_der=.1;end;
frwx=logspace(log10(5.e-4),log10(300),500);


u2y_sys=input('Give system file  ','s');

sch_tol=input('Give reduction tolerance  [1.e-6]');
if length(sch_tol)<1;sch_tol=1.e-6;end

eval(['load ' u2y_sys])
SI_=mksys(Apr,Bpr,Cpr,Dpr,'ss');
[Adi,Bdi,Cdi,Ddi]=ssinv(Fpr,Bypr,-Cpr,eye(size(Dypr))-Dypr);
D_P=mksys(Adi,Bdi,Cdi,Ddi,'ss');
N_P=mksys(Fpr,Bpr-Bypr*Dpr,Cpr,Dpr,'ss');

[Noutp,Ninp]=size(Dpr);
if Noutp~=Ninp
   disp('pid_des cannot handle non-square systems')
end

disp('Step 1:  extract inner diagonals and compute Perron sv')
per_tem=input('Do you want to skip this step? [0=no]  ');
if isempty(per_tem);per_tem=0;end
if per_tem ~=1
   [Apn,Bpn,Cpn,Dpn,Noutp,Ninp]=extr_dia(Apr,Bpr,Cpr,Dpr,sch_tol);
   [Ax,Bx,Cx,Dx]=addss(Apr,Bpr,Cpr,Dpr,Apn,Bpn,-Cpn,-Dpn);
   svPx=sv3_5(Ax,Bx,Cx,Dx,-1,frun);
   svPxs=sv3_5(Ax,Bx,Cx,Dx,1,frun);svPxs=svPxs(1,:);
   svPx=min([svPx;svPxs]);
   svPP=sv3_5(Apn,Bpn,Cpn,Dpn,1,frun);
   MUNCx=(svPx)'./min(svPP)';
   clf
   loglog(frun,sat(MUNCx,1.e4),frun,sat(MUNC(:,1)+MUNCx,1.e4),'--');grid;
   title('Approximation Error SV'); pause
end
%---------------------------------------------------------
disp('Step 2:  TUNE PIDs ')
tc_filti=zeros(Noutp,1);KPIDIN=zeros(3*Noutp,3);tdone=0;
while tdone==0
   i_mult=ones(1,Ninp).*i_mult1;i_mult=diag(i_mult);
   Bpr=Bpr*i_mult;Dpr=Dpr*i_mult;
   lwdef=input('Give default value for freq.weight [1]  ');
   if length(lwdef)<1;lwdef=1;end
   
   for IX = 1:Noutp
      disp('Tuning loop:'),disp(IX)
      rdone=input('Done with this loop? [0=n]  ');
      if isempty(rdone);rdone=0;end
      if rdone == 0
         %------  MODEL REDUCTION --------
         A_sys=mksys(Apr,Bpr(:,IX),Cpr(IX,:),Dpr(IX,IX),'ss');
         [Aprx,Bprx,Cprx,Dprx]=w_sysred(A_sys,[],[],0,[],[3,sch_tol]);
         %------  TUNING -----------------
         disp('roots of numerator and denominator')
         [NUMx,DENx]=ss2tf(Aprx,Bprx,Cprx,Dprx,1);
         roots(NUMx),roots(DENx)
         T_fil=input('Error Filter TC [0.1]   ');
         if length(T_fil)<1;T_fil=.1;end;pidfil=[T_fil 1];
         tc_filti(IX)=T_fil;
         
         disp('Tuner options: type 0 = 1st order loopshape (default)')
         disp('               type 1 = 2nd order loopshape (integrator)')
         disp('               type 2 = 2nd order loopshape (pole-zero auto select)')
         disp('               type 3 = 2nd order loopshape (pole-zero user defined)')
         disp('Recom: Use optimization method = 1.1 with type 1')
         disp('Format: [BW,(type,zero,pole)]; auto: pole=slowest, zero=zero*BW')
         
         lam1=input('give bandwidth and type of PID design  ');
         op_t=input('give optimization method 0=Hinf, val>1=H2/Hinf [0]  ');
         if length(op_t)<1;op_t=0;end
         if length(lam1)<1;lam1=1;end
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
   end
   %--------  EVALUATING DESIGN  ------------------
   Bpr=Bpr*inv(i_mult);Dpr=Dpr*inv(i_mult);
   
   kpidin=KPIDIN([3:3:3*Noutp],:);
   kpidin(:,1)=i_mult*kpidin(:,1);
   
   [Acr,Bcr,Ccr,Dcr]=npidss(kpidin,T_der,tc_filti,0);
   C_I=mksys(Acr,Bcr,Ccr,Dcr,'ss');
   [Algn,Blgn,Clgn,Dlgn]=series(Acr,Bcr,Ccr,Dcr,Apr,Bpr,Cpr,Dpr);
   [Acly,Bcly,Ccly,Dcly]=feedbk(Algn,Blgn,Clgn,Dlgn,2);
   clf;TIM=[0:T_max/300:T_max]';YIS=zeros(length(TIM),Noutp*Noutp);
   for IX=1:Noutp
      YIS(:,(IX-1)*Noutp+1:IX*Noutp)=step(Acly,Bcly,Ccly,Dcly,IX,TIM);
   end
   plot(TIM,YIS);title('Step Responses');pause
   Y=lsim(Acly,Bcly,Ccly,Dcly,ones(length(TIM),Ninp),TIM);plot(TIM,Y);
   title('All-channel Step Response');pause
   [LPMG] = abs(freqrc(Acly,Bcly,Ccly,Dcly,frun));
   loglog(frun,LPMG);title('Frequency Responses');pause
   %-------------------------------------------- robustness checks
   disp('** Checking Robust stability of the spike loop design...')
   per_tem=input('Do you want to skip this step? [0=no]  ');
   if isempty(per_tem);per_tem=0;end
   if per_tem ~=1
      [stchk,MUNCe,S_cl]=unc_chkg(N_P,D_P,C_I,frun,AUNC,flag);
      if max(stchk) > .95
         disp('*****  Potential robustness problems in the loop!!')
         flag=max(flag,1);
      else
         disp('*****  Design seems successful!!')
      end
      pause;loglog(frun,stchk);grid;title('Robust stability condition');pause
   end
   %---------------------------------------------------------------
   tdone=input('Done?  [1=yes]  ');
   if isempty(tdone);tdone=1;end
end
savsys=input('Give filename to save results (e.g., out...)  ','s');
if length(savsys)>2
   eval(['save ',savsys,' kpidin tc_filti T_der Acr Bcr Ccr Dcr'])
end
%-------------------------------------------- robustness checks

