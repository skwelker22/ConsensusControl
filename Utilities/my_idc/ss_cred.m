function [Acr,Bcr,Ccr,Dcr,scc]=s_cred(AC_,AP_,AW_,frw,flag,cl_spec);
%function [Acr,Bcr,Ccr,Dcr]=s_cred(AC_,AP_,AW_,frw,flag,cl_spec);
% Compensator reduction.
%    Slow-fast decompositions are used to remove too fast
%    and too slow dynamics; 
%    defaults are: 100 x T_BW (cl_spec(1)),
%                  0.03 x S_BW (cl_spec(3)).
%  AC_ is the packed full order compensator
%  AP_ is the plant 
%  AW_ contains plant augmentation (e.g., integrators)
%  frw is a frequency vector for the computation and display 
%     of singular values; 0 computes SV's 2 dec above and below T_BW;
%     displays are active for flag>1.

scc=1;
if nargin < 5;flag=1;end
if nargin < 6;cl_spec=1;end
if length(cl_spec)<3;
   cl_spec=[cl_spec(1),0,cl_spec(1)/2,0,1,1];
end
CUT_HI=cl_spec(1)*100;CUT_LO=cl_spec(3)/30;
if frw==0;
   frw=logspace(log10(cl_spec(1)/100),log10(cl_spec(1)*100),100)';
end
[acp,bcp,ccp,dcp]=branch(AC_);
[Apr,Bpr,Cpr,Dpr]=branch(AP_);
[Awa,Bwa,Cwa,Dwa]=branch(AW_);
[Apa,Bpa,Cpa,Dpa]=series(Awa,Bwa,Cwa,Dwa,Apr,Bpr,Cpr,Dpr);
EAC=eig(acp);
if flag>=2
   magu=sv3_5(acp,bcp,ccp,dcp,1,frw);
   loglog(frw,magu);title('Compensator SV (without augmentation)')
   if flag<=4;pause(2);else;pause;end
end
if flag>3
   disp('Magnitude of Compensator Poles')
   abs(EAC)'
   disp('Magnitude of Compensator Transmission Zeros')
   TZAC=tzero(acp,bcp,ccp,dcp);abs(TZAC)'
   disp('Maximum real parts of poles and zeros')
   disp([max(real(EAC)),max(real(TZAC))])
end
%  -------SLOW-FAST DECOMPOSITION (HI-FREQ)------------
if flag > 1
  CUTfr=input(['cutoff freq  [',num2str(CUT_HI),']   ']);
else
  CUTfr=[];
end

if length(CUTfr)<1;CUTfr=CUT_HI;end
CUT=length(find(abs(EAC)<=CUTfr))
if CUT == 0; CUT=length(acp);end
  if CUT < length(acp)
    [Acr,Bcr,Ccr,Dcr,AH,BH,CH,DH] = slowfast(acp,bcp,ccp,dcp,CUT);
    disp('Including hf-DC-component in D')
    DH-CH*inv(AH)*BH;
        Dcr=Dcr+DH-CH*inv(AH)*BH;
  else
    Acr=acp;Bcr=bcp;Ccr=ccp;Dcr=dcp;
  end

   if flag>3
     magur=sv3_5(Acr,Bcr,Ccr,Dcr,1,frw);
     loglog(frw,magu,frw,magur,'--')
     title('Controller singular values'),
     if flag<=4;pause(2);else;pause;end
   end
     scc=h_perfev(Acr,Bcr,Ccr,Dcr,Apa,Bpa,Cpa,Dpa,frw,flag-.5);
     
%  -------SLOW-FAST DECOMPOSITION (LOW-FREQ)------------
if scc ~=0
   disp('Low freq. reduction')
     if flag >2
       itemp=input('Reject low frequencies (0=no)? [1] ');
       if isempty(itemp);itemp=1;end
     else
       itemp=1;
     end
else
   itemp=0;
end

if itemp ~= 0
   EAC=eig(Acr);
     if flag>1
        disp('Magnitude of Compensator Poles')
        abs(EAC)'
     end
     if flag >1
        CUTfr=input(['cutoff freq  [',num2str(CUT_LO),']   ']);
     else
        CUTfr=[];
     end
   if length(CUTfr)<1;CUTfr=CUT_LO;end
   CUT=length(find(abs(EAC)<=CUTfr))
   if CUT==0;CUT=length(Acr);end
     if CUT < length(Acr)
        [AH,BH,CH,DH,Acr,Bcr,Ccr,Dcr] = slowfast(Acr,Bcr,Ccr,Dcr,CUT);
        disp('Including hf-component in D')
        Dcr=Dcr+DH;
     end
   if flag>3
     magur=sv3_5(Acr,Bcr,Ccr,Dcr,1,frw);
     loglog(frw,magu,frw,magur,'--')
     title('Controller singular values'),
     if flag<=4;pause(2);else;pause;end
   end
   scc=h_perfev(Acr,Bcr,Ccr,Dcr,Apa,Bpa,Cpa,Dpa,frw,flag-.5);
end

%---------------------WEIGHTED CONTROLLER REDUCTION------
if scc ~= 0
   if flag >2
      itemp=input('Perform Weighted Reduction? (0=n) [1]  ');
      if isempty(itemp);itemp=1;end
   else
      itemp=1;
   end
else
   itemp=0;
end

if itemp~=0
  done=0;flagt=flag;
  [Aw,Bw,Cw,Dw]=feedbk(Apa,Bpa,Cpa,Dpa,3,Acr,Bcr,Ccr,Dcr);
  W_=mksys(Aw,Bw,Cw,Dw,'ss');C_=mksys(Acr,Bcr,Ccr,Dcr,'ss');

    while done==0
          if flagt > 1
            [Acr,Bcr,Ccr,Dcr]=w_sysred(C_,[],W_,flagt,0,0);
          else
            [Acr,Bcr,Ccr,Dcr]=w_sysred(C_,[],W_,flagt,0,[3,5.e-5]);
          end
      scc=h_perfev(Acr,Bcr,Ccr,Dcr,Apa,Bpa,Cpa,Dpa,frw,flag+2);
          if flagt >=1
            done=input('Done with weighted reduction? (0=no) [1]  ');
            if isempty(done);done=1;end
            if done==0;flagt=2;end
          else
            done=1;
          end
    end
end


