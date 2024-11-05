function [Acd,Bcd,Ccd,Dcd]=m_d_conu(P_u,C_y,P_d,W_d,lam,G_w,e_s,flag);
% function [Acd,Bcd,Ccd,Dcd]=m_d_conu(P_u,C_y,P_d,W_d,lam,G_w,e_s,flag);
% feedforward controller for a measured disturbance (input matching)
%   P_u, C_y: Packed systems, u->y plant, y->u controller
%   P_d: dist -> y (packed)
%   W_d: biproper disturbance weight (e.g., frequency content)
%        stable and min-phase.
%   lam: scalar weight penalizing corrections
%        makes the problem well defined for more ins than outs.
%        1e-4 if not specified (empty).
%   G_w: packed bias term (ninp*ndis; 0 if not specified)
%   e_s: [exp_shift, low-fr cutoff, Nehari=0/H2=1]
% On return
%   [Acd,Bcd,Ccd,Dcd] or
%   C_dis: packed controller; u = C_y (r-y) - C_dis (d)

% KST 9/98

frw=logspace(-4,4,100)';

if nargin<4;W_d=[];end
if nargin<5;lam=[];end;  if isempty(lam);lam=1e-4;end
if nargin<6;G_w=[];end
if nargin<7;e_s=[1e-5,1e-3,0];end; 
if isempty(e_s);e_s=[1e-5,1e-3,0];end
if length(e_s)<2;e_s=[e_s,1e-3];end
if length(e_s)<3;e_s=[e_s,0];end
if nargin<8;flag=1;end
[Ap,Bp,Cp,Dp]=branch(P_u);
[Ac,Bc,Cc,Dc]=branch(C_y);
[Ad,Bd,Cd,Dd]=branch(P_d);
[noutp,ninp]=size(Dp);
[ntem,ndis]=size(Dd);
if ~isempty(W_d);[Aw,Bw,Cw,Dw]=branch(W_d);end
if ~isempty(G_w);[Ag,Bg,Cg,Dg]=branch(G_w);end

% create sensitivities
[Alp,Blp,Clp,Dlp]=series(Ap,Bp,Cp,Dp,Ac,Bc,Cc,Dc);
[At2,Bt2,Ct2,Dt2]=feedbk(Alp,Blp,Clp,Dlp,2);
[Asp,Bsp,Csp,Dsp]=feedbk(Ac,Bc,Cc,Dc,3,Ap,Bp,Cp,Dp);
[At1,Bt1,Ct1,Dt1]=series(Ad,Bd,Cd,Dd,Asp,Bsp,Csp,Dsp);

if flag>=2
   sr1=sv3_5(At2,Bt2,Ct2,Dt2,1,frw);
   sr2=sv3_5(At2,Bt2,Ct2,Dt2+lam*eye(size(Dt2)),1,frw);
   loglog(frw,sr1,'r',frw,sr2,'g');pause
end

% create nehari tfms
disp('--- converting to Nehari...');
Iinp=eye(ninp,ninp);Idis=eye(ndis,ndis);
if lam>0
   [U,S,V]=svd(Dt2);
   SD1=diag(max(diag(S),lam));nsd1=length(SD1);
   S(1:nsd1,1:nsd1)=SD1;
   Dt2=U*S*V';
   [U,S,V]=svd(Dt1);
   SD1=diag(max(diag(S),lam));nsd1=length(SD1);
   S(1:nsd1,1:nsd1)=SD1;
   Dt1=U*S*V';
end
if ~isempty(W_d);
   [At1,Bt1,Ct1,Dt1]=series(Aw,Bw,Cw,Dw,At1,Bt1,Ct1,Dt1);
end

% shift slightly to avoid zeros on jw axis
At1s=At1-e_s(1)*eye(size(At1));
At2s=At2-e_s(1)*eye(size(At2));
if ~isempty(W_d);Aws=Aw-e_s(1)*eye(size(Aw));end

% done with setup:  min_Cd ||T1 - T2*Cd*W||
% convert to nehari min_X ||R - X ||;     X = T2o*Cd*W
disp('--- Computing IOFR...');

[AI,BI,CI,DI,AP,BP,CP,DP,AO,BO,CO,DO]=iofr(At2s,Bt2,Ct2,Dt2);
AI_=mksys(AI,BI,CI,DI,'ss');
%AO_=mksys(AO,BO,CO,DO,'ss');
[AI,BI,CI,DI]=h_sysred(AI_,[],[],flag-2,[],[3,1e-4]);
%[AO,BO,CO,DO]=h_sysred(AO_,[],[],flag,[],0);
%  AIa=-AI';BIa=-CI';CIa=BI';DIa=DI';
[AR,BR,CR,DR]=series(At1s,Bt1,Ct1,Dt1,-AI',-CI',BI',DI');

disp('--- Solving Nehari...');
if e_s(3)==0
   [AX,BX,CX,DX,sigr]=nehari(AR,BR,CR,DR);
   disp(['   Nehari min distance = ',num2str(sigr)]);  
else
   [AX,BX,CX,DX,xx1,xx2,xx3,xx4] = stabproj(AR,BR,CR,DR);
   DX=DX+xx4;
end
if flag>=2
%sr=sv3_5(AR,BR,CR,DR,1,frw);
%sx=sv3_5(AX,BX,CX,DX,1,frw);
   [Arx,Brx,Crx,Drx]=addss(AR,BR,CR,DR,AX,-BX,CX,-DX);
   sx=sv3_5(Arx,Brx,Crx,Drx,1,frw);
   loglog(frw,sx);pause
end
%AX=AX+e_shift*eye(size(AX));

% Solve for Cd by inverting the outer factors of X
%   W is square and invertible

%if ~isempty(W_d);
%   [Awi,Bwi,Cwi,Dwi]=ssinv(Aws,Bw,Cw,Dw);
%   [Acd,Bcd,Ccd,Dcd]=series(Awi,Bwi,Cwi,Dwi,AX,BX,CX,DX);
%else
%   Acd=AX;Bcd=BX;Ccd=CX;Dcd=DX;
%end

Acd=AX;Bcd=BX;Ccd=CX;Dcd=DX;

DOi=pinv(DO);
AOi=AO-BO*DOi*CO;
BOi=BO*DOi;COi=-DOi*CO;
[Acd,Bcd,Ccd,Dcd]=series(Acd,Bcd,Ccd,Dcd,AOi,BOi,COi,DOi);
Acd=Acd+e_s(1)*eye(size(Acd));

disp('--- Performing reductions...');
%  -------SLOW-FAST DECOMPOSITION (LOW-FREQ)------------
if flag >2
   itemp=input('Reject low frequencies (0=no)? [1] ');
   if isempty(itemp);itemp=1;end
else
   itemp=1;
end
if itemp ~= 0
   EAC=eig(Acd);
   if flag>1
      disp('Magnitude of Compensator Poles')
      abs(EAC)'
   end
   if flag >1
      CUTfr=input(['cutoff freq  [',num2str(e_s(2)),']   ']);
   else
      CUTfr=[];
   end
   if length(CUTfr)<1;CUTfr=e_s(2);end
   CUT=length(find(abs(EAC)<=CUTfr))
   if CUT==0;CUT=length(Acd);end
   if CUT < length(Acd)
      [AH,BH,CH,DH,Acd,Bcd,Ccd,Dcd] = slowfast(Acd,Bcd,Ccd,Dcd,CUT);
      disp('Including hf-component in D')
      Dcd=Dcd+DH;
   end
end

C_dis=mksys(Acd,Bcd,Ccd,Dcd,'ss');
[Acd,Bcd,Ccd,Dcd]=h_sysred(C_dis,W_d,[],flag-2,[],[3,1e-4]);

%if ~isempty(W_d)
%   [Acx,Bcx,Ccx,Dcx]=series(Aw,Bw,Cw,Dw,Acd,Bcd,Ccd,Dcd);
%else
%   Acx=Acd;Bcx=Bcd;Ccx=Ccd;Dcx=Dcd;
%end

Acx=Acd;Bcx=Bcd;Ccx=Ccd;Dcx=Dcd;
[Acx,Bcx,Ccx,Dcx]=series(Acx,Bcx,Ccx,Dcx,At2,Bt2,Ct2,Dt2);

if flag>=2
   sr=sv3_5(At1,Bt1,Ct1,Dt1,1,frw);
   sx=sv3_5(Acx,Bcx,Ccx,Dcx,1,frw);
   loglog(frw,sr,'r',frw,sx,'b');pause
end

if nargout==1
   C_dis=mksys(Acd,Bcd,Ccd,Dcd,'ss');
   Acd=C_dis;
end