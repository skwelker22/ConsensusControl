function [Acr,Bcr,Ccr,Dcr,savsys,cl_spco]=s_hinf(Apr,Bpr,Cpr,Dpr,frun,MUNC,flag,cl_spec,fl_ns,figno);

%function [Acr,Bcr,Ccr,Dcr,savsys,cl_spco]=
%             s_hinf(Apr,Bpr,Cpr,Dpr,[frun,MUNC,flag,cl_spec,fl_ns,figno]);
%         for the design of an H-inf controller.
%   The plant should be in [Apr,Bpr,Cpr,Dpr]
%   Diagonal performance/robustness weights are specified through 
%   cl_spec, or interactively; 
%   Uncertainty estimates should be in [frun,MUNC] to guide the weight
%     selection; otherwise use [] 
%   An exponential shift is used for the plant only
%     so that no poles are too close to the jw-axis (flag>1); 
%   Diagonal plant augmentation by integrators [s+a]/s can
%     be specified (flag>1) or 1/s is used as default

format short e
disp('------------------------------  Function S_hinf executing....')


if nargin < 5;frun=[];MUNC=[];flag=0;end
if nargin < 7; flag=0;end
if nargin < 8; cl_spec=[];end
if nargin < 9; fl_ns=[];end
if nargin < 10; figno=1;end

if length(cl_spec)>=1;cl_bw=cl_spec(1);else;cl_bw=3;end
frw=logspace(log10(cl_bw/100),log10(cl_bw*100),100)';
if length(frun)==0;frun=frw;end
if length(MUNC)==0;MUNC=[frun,frun]*0+1;end

%----------------------------------------------------------------
[noutp,ninp]=size(Dpr);
%----------------------------------------------------------------
%----------------------------------- Integrator Augmentation ----
if flag>1
  disp('Plant Eigenvalues')
  eig(Apr)
  disp('Plant transmission zeros')
  tzero(Apr,Bpr,Cpr,Dpr)
end
if flag >1
  itemp=input('Augment with integrators? (1,2,...=yes) [1]  ');
else
  itemp=[];
end
numaw=[];
if length(itemp) >0;int_no=itemp;else;int_no=1;end
if int_no>0 & flag >=1
   numaw=input(['Give Numerator Polynomial [1 ',num2str(cl_bw),']  ']);
   if length(numaw)<1;numaw=[1 cl_bw];end
elseif int_no>0 & flag <=1
   numaw=[1,cl_bw];
end
if length(numaw)>0 
     [Awa,Bwa,Cwa,Dwa]=hinf_wgt(ninp,numaw,[1 zeros(1,int_no)],0);
else
   Awa=[];Bwa=[];Cwa=[];Dwa=eye(ninp,ninp);
end

       AW_=mksys(Awa,Bwa,Cwa,Dwa,'ss');
       AP_=mksys(Apr,Bpr,Cpr,Dpr,'ss');
[Apra,Bpra,Cpra,Dpra]=series(Awa,Bwa,Cwa,Dwa,Apr,Bpr,Cpr,Dpr);
Dprad=Dpra;
%----------------------------------------------------------------
%---------------------------------------- Weight Specs ----------
   if flag>=1
     disp('S_hinf specs [T_bw,T_roll,S_bw,S_roll,RISK,AGRS]')
     disp(cl_spec)
     bwdth=input('s_hinf: Specify target loop T,S ? (1=yes) [0]  ');
     if isempty(bwdth);bwdth=0;end
     if length(cl_spec)<4;bwdth=1;disp('YES!!!');end
     if bwdth==1
       disp('Enter Robustness Weight');
       [Po3,Ze3,Ga3,ftem3,mtem3]=wgt_sel(noutp,frun,MUNC,0,flag);
       cl_spec(1)=min(Po3);cl_spec(2)=length(Po3);
       cl_bw=cl_spec(1);
       [Aw3,Bw3,Cw3,Dw3]=hinf_wgt(noutp,Po3,Ze3,Ga3);
       disp('Enter Performance Weight');
       [Po1,Ze1,Ga1,ftem1,mtem1]=...
           wgt_sel(noutp,ftem3,[(1)./mtem3 MUNC],1,flag);
       cl_spec(3)=max(Po1);cl_spec(4)=length(Po1);
       [Aw1,Bw1,Cw1,Dw1]=hinf_wgt(noutp,Po1,Ze1,Ga1);
     else
       ftem3=frun;ftem1=frun;
       [mtem3,num3,den3]=awgt_sel(ftem3,cl_spec(1),cl_spec(2),0);
       [mtem1,num1,den1]=awgt_sel(ftem1,cl_spec(3),cl_spec(4),-1);
       [Aw3,Bw3,Cw3,Dw3]=hinf_wgt(noutp,den3,num3,0);
       [Aw1,Bw1,Cw1,Dw1]=hinf_wgt(noutp,den1,num1,0);
     end
   else
       if length(cl_spec)<6
          cl_spec=q_spec(frun,MUNC);
       end
       ftem3=frun;ftem1=frun;
       [mtem3,num3,den3]=awgt_sel(ftem3,cl_spec(1),cl_spec(2),0);
       [mtem1,num1,den1]=awgt_sel(ftem1,cl_spec(3),cl_spec(4),-1);
       [Aw3,Bw3,Cw3,Dw3]=hinf_wgt(noutp,den3,num3,0);
       [Aw1,Bw1,Cw1,Dw1]=hinf_wgt(noutp,den1,num1,0);
   end

  MTMP2=sat(mtem1,100,1.e-5);MTMP3=sat(mtem3,100,1.e-5);
  MTMP4=sat(1.0./MUNC,100,1.e-5);
figure(figno);clf
  loglog(ftem3,MTMP3,ftem1,MTMP2,frun,MTMP4);grid
  title('Weights and preliminary uncertainties')

Aw2=[];Bw2=[];Cw2=[];Dw2=eye(ninp,ninp);
%-------------------------------------------------------------------
%---------------------------------- H-inf solution parameters ------
   if flag>=1
      disp('give: 1. plant exp.shift, 2. gain of contr.weight, ')
      disp('      3. perf_only switch, in [ ], or hit return for defaults ')
      disp('for example: [(.001*cl_bw),.01,0] ')
      spectemp=input(' ');
   else
      spectemp=[];
   end

cgain=.01; dlev=.001*cl_bw;perfonly=0;
dshift=0;pgain=1; rgain=1;gamaux=[.005 1.2 0];
if length(spectemp)>0,dlev=spectemp(1);end
if length(spectemp)>1,cgain=spectemp(2);end
if length(spectemp)>2,perfonly=spectemp(3);end
spectemp=[dlev,cgain,perfonly];
   if flag>1;disp(spectemp);end

LP_BW=sv3_5(Apra,Bpra,Cpra,Dpra,1,cl_bw);LP_BW=min(LP_BW);
thrbi=1.e-5*abs(LP_BW);
     PLSHP0=sv3_5(Apra,Bpra,Cpra,Dpra,1,frw);
     PLSWT=sv3_5(Aw3,Bw3,Cw3,Dw3,1,frw);PLSWT=1./PLSWT;
if flag>=1;itemp=0;else;itemp=1;end
[UD,SD,VD]=svd(Dprad);
SD=diag(diag(SD)+thrbi);Dpra=UD*SD*VD';
diag(SD)'
   while itemp ==0
figure(figno+1);clf
      PLSHP=sv3_5(Apra,Bpra,Cpra,Dpra,1,frw);
      loglog(frw,PLSHP,frw,PLSHP0,frw,PLSWT,'g:');
      title('Augmented plant singular values and weight');
      itemp=input('Done with thruput fix? (0=no) [1]  ');
      if isempty(itemp);itemp=1;end
         if itemp==0
           thrbi=input('biproperness threshold [1.e-4]  ');
           if length(thrbi)<1,thrbi=1.e-4;end
           [UD,SD,VD]=svd(Dprad);
           SD=diag(diag(SD)+thrbi);Dpra=UD*SD*VD';
         end
   end
%---------------------------------------------------------------
%----------------------- augm. plant Apra,Bpra,Cpra,Dpra -------
%-------------------------------- Super-plant construction -----
nsp=length(Apra);Apra=Apra+dlev*eye(nsp,nsp);
Cw1=Cw1*pgain;Dw1=Dw1*pgain;
Cw2=Cw2*cgain;Dw2=Dw2*cgain;
Cw3=Cw3*rgain;Dw3=Dw3*rgain;
ns1=length(Aw1);ns2=length(Aw2);ns3=length(Aw3);
A=[Apra, zeros(nsp,ns1+ns2+ns3)];
A=[A;-Bw1*Cpra,Aw1,zeros(ns1,ns2+ns3)];
A=[A; zeros(ns2,nsp+ns1),Aw2,zeros(ns2,ns3)];
A=[A; Bw3*Cpra,zeros(ns3,ns1+ns2),Aw3];
A=A+dshift*eye(nsp+ns1+ns2+ns3,nsp+ns1+ns2+ns3);
B1=[zeros(nsp,noutp);Bw1;zeros(ns2,noutp);zeros(ns3,noutp)];
B2=[Bpra;-Bw1*Dpra;Bw2;Bw3*Dpra];
C1=[-Dw1*Cpra,Cw1,zeros(noutp,ns3+ns2)];
C1=[C1;zeros(ninp,nsp+ns1),Cw2,zeros(ninp,ns3)];
C1=[C1;Dw3*Cpra,zeros(noutp,ns1+ns2),Cw3];
C2=[-Cpra,zeros(noutp,ns1+ns3+ns2)];
D11=[Dw1;zeros(ninp,noutp);zeros(noutp,noutp)];
D12=[-Dw1*Dpra;Dw2;Dw3*Dpra];
D21=[eye(noutp,noutp)];D22=-Dpra;
[ny1,nu1]=size(D11);[ny2,nu2]=size(D22);
BB=[B1,B2];CC=[C1;C2];DD=[D11,D12;D21,D22];
%----------------------------------------------------
eimax=max(real(eig(A)));
if flag > 1
  disp(['Max real eigenvalue of super-plant = ',num2str(eimax)])
  itemp=0; disp('-----Model Order Reduction-----')
    while itemp==0
      delred=input('give delta-shift for mod-ord red. []  ')
         if length(delred)<1;delred=max(eimax+0.015*cl_bw,0);end
      Atem=A-delred*eye(length(A),length(A));
      hsv=hksv(Atem,BB,CC);
      disp('Super-plant Hankel singular values'),hsv'
      itemp=input('DONE with shifting? (0=no) [1] ');
      if length(itemp)<1;itemp=1;end
    end
  nord=input('give order of reduced model [all] ');
  if isempty(nord);nord=length(Atem);end
    if nord < length(hsv)
      [AM BM CM DM totbnd hsv]=schmr(Atem,BB,CC,DD,1,nord);
    else 
      AM=Atem;BM=BB;CM=CC;DM=DD;
    end
else
  delred=max(eimax+0.02*cl_bw,0);
  Atem=A-delred*eye(length(A),length(A));
  AM=Atem;BM=BB;CM=CC;DM=DD;
end
% -------------------------- Super-plant Balancing ----------
   if flag >1
      ibalan=input('Perform balancing (0=no) [1] ');
      if isempty(ibalan);ibalan=1;end
   else
     ibalan=1;
   end
   if ibalan ~=0
     [AM,BM,CM,tem1,tem2]=obalreal(AM,BM,CM);
   end
A=AM+delred*eye(length(AM),length(AM));
B1=BM(:,1:nu1);B2=BM(:,nu1+1:nu1+nu2);
C1=CM(1:ny1,:);C2=CM(ny1+1:ny1+ny2,:);
D11=DM(1:ny1,1:nu1);D12=DM(1:ny1,nu1+1:nu1+nu2);
D21=DM(ny1+1:ny1+ny2,1:nu1);
D22=DM(ny1+1:ny1+ny2,nu1+1:nu1+nu2);
% ------------------------------- pack matrices ------------------
TSS_ = mksys(A,B1,B2,C1,C2,D11,D12,D21,D22,'tss');
  if perfonly==1
    gamindex=0*[1:ny1];
    gamindex(1:noutp)=gamindex(1:noutp)+1;
  else
    gamindex=[];
  end
%----------------------------------- H-inf solution ----------

[gamaopt,SS_CP]=hinfopt(TSS_,gamindex,gamaux);

% --------------------------------- unpack controller --------
[acp,bcp,ccp,dcp]=branch(SS_CP);
acp=acp-(dshift+dlev)*eye(length(acp),length(acp));
   if flag >=2
     disp('Controller eigenvalues '),eig(acp)
   end
% -------------------------------- Compute Sensitivities -------
[Alpgn,Blpgn,Clpgn,Dlpgn]=series(acp,bcp,ccp,dcp,...
            Apra-(dshift+dlev)*eye(nsp,nsp),Bpra,Cpra,Dpra);
[Acly,Bcly,Ccly,Dcly]=feedbk(Alpgn,Blpgn,Clpgn,Dlpgn,2);
[Acle,Bcle,Ccle,Dcle]=feedbk(Alpgn,Blpgn,Clpgn,Dlpgn,1);
magcl=sv3_5(Acly,Bcly,Ccly,Dcly,1,frw);
magcle=sv3_5(Acle,Bcle,Ccle,Dcle,1,frw);
%magw1=sv3_5(Aw1,Bw1,Cw1,Dw1,1,frw);magw3=sv3_5(Aw3,Bw3,Cw3,Dw3,1,frw);
MTMP0=sat(magcl',100,1.e-5);MTMP1=sat(magcle',100,1.e-5);
MTMP2=sat(mtem1,100,1.e-5);MTMP3=sat(mtem3,100,1.e-5);
MTMP4=sat(1.0./MUNC,100,1.e-5);

figure(figno+2);clf
loglog(frw,MTMP0,'g',frw,MTMP1,'b',ftem1,MTMP2,'y',ftem3,MTMP3,'y',frun,MTMP4,'r')
title('Cl.Lp. T (g), S (b) Singular Values; (weights=y, UNC=r)')
max_T=max(max(magcl));max_S=max(max(magcle));
disp(['max |T| = ',num2str(max_T),',   max |S| = ',num2str(max_S)])
   if max(real(eig(Acly)))<0;scc=1;else;scc=0;end

if flag >=1
   itemp=input('Reduce compensator? (0=no) [1] ');
   if isempty(itemp);itemp=1;end
else
  itemp=1;
end
   if itemp ~=0;rdone=0;flagt=flag;
      while rdone==0
         [Acr,Bcr,Ccr,Dcr,scc]=s_cred(SS_CP,AP_,AW_,frw,flagt,cl_spec,figno+3);
            if flagt>=1
               rdone=input('Done with reduction? (0=no) [1]  ');
               if isempty(rdone);rdone=1;end
            else
              if scc ~=0;rdone=1;else;rdone=0;end
            end
            
         if rdone==0;flagt=10;end
      end
   else
      [Acr,Bcr,Ccr,Dcr]=series(acp,bcp,ccp,dcp,Awa,Bwa,Cwa,Dwa);
      scc=0;
      fl_ns =[];
   end

R_A=inv_spec(cl_spec,frun,MUNC);
if length(cl_spec)<6
   cl_spec(5:6)=R_A;
else
   if abs(log(R_A(1)/cl_spec(5)))+abs(R_A(2)-cl_spec(6))>0.002
     cl_spec(5:6)=R_A;
   end
end
cl_spco=cl_spec;

if scc ~=0
   disp('Nominal H-inf controller design successful!')
else
   disp('NOMINAL H-INF CONTROLLER DESIGN UNSUCCESSFUL!')
   fl_ns=[];
end

if length(fl_ns)>=2 & flag <=1 
   savsys=fl_ns;
else
   savsys=input('Give filename to save results (e.g., out...)  ','s');
end
if length(savsys)>=2
  eval(['save ',savsys,' SS_CP AP_ AW_ cl_spec scc Acr Bcr Ccr Dcr'])
end

