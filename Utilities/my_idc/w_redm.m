function [Acr,Bcr,Ccr,Dcr]=w_redm(AC_,AWi_,AWo_,def_flag,thr,nn);
%function [Acr,Bcr,Ccr,Dcr]=w_redm(AC_,AWi_,AWo_,def_flag,thr,nn);
% Compensator reduction by weighted balanced truncation method.
%     "min" || Wo (C-Cr) Wi ||_inf
%  AC_ is the packed full order compensator (assumed minimal)
%  AWi_ is the input weight and AWo_ is the output weight
%         (stable)
%  thr = [i o] biproperness thresholds for i,o weights
%  nn  = [type,n]: reduction type and cutoff point
%         type = 0: display SV and ask for n
%         type = 1: keep the n highest SV
%         type = 2: keep SV's higher than n
%         type = 3: keep SV's larger than n*maxSV

if nargin < 4;def_flag=1;end
if nargin < 5;thr=0;end
if isempty(thr);thr=0;end
if thr==0;thr=[1 1]*1.e-4;end
if length(thr)==1;thr=[1 1]*thr;end
if nargin<6;nn=0;end
if nn(1)~=0 & length(nn)<2; nn=0;end

[Ac,Bc,Cc,Dc]=branch(AC_);
if max(real(eig(Ac)))>=0
   [AC_s,AC_a] = stabproj(AC_);
   [Aca,Bca,Cca,Dca]=branch(AC_a);
   [Ac,Bc,Cc,Dc]=branch(AC_s);nc=length(Ac);
else
   nc=length(Ac);
   Aca=[];Bca=[];Cca=[];Dca=0*Dc;
end

if AWi_==[]
   iflag=1;Ai=[];Bi=[];Ci=[];Di=eye(size(Dc'*Dc));
else
   iflag=0;[Ai,Bi,Ci,Di]=branch(AWi_);
end
if AWo_==[]
   oflag=1;Ao=[];Bo=[];Co=[];Do=eye(size(Dc*Dc'));
else
   oflag=0;[Ao,Bo,Co,Do]=branch(AWo_);
end
Dc0=Dc;Dc=0*Dc;
         [UD,SD,VD]=svd(Di);
         SD1=diag(max(diag(SD),thr(1)));nsd1=length(SD1);
         SD(1:nsd1,1:nsd1)=SD1;
         Di=UD*SD*VD';
         [UD,SD,VD]=svd(Do);
         SD1=diag(max(diag(SD),thr(2)));nsd1=length(SD1);
         SD(1:nsd1,1:nsd1)=SD1;
         Do=UD*SD*VD';
if iflag==0
  [Ani,Bni,Cni,Dni,Ax,Bx,Cx,Dx,Ai,Bi,Ci,Di] = iofc(Ai,Bi,Ci,Di); 
end
if oflag==0
  [Ami,Bmi,Cmi,Dmi,Ax,Bx,Cx,Dx,Ao,Bo,Co,Do] = iofr(Ao,Bo,Co,Do); 
end

[Ax,Bx,Cx,Dx]=series(-Ai,Bi,-Ci,Di,Ac,Bc,Cc,Dc);
[Ax,Bx,Cx,Dx]=series(Ax,Bx,Cx,Dx,-Ao,Bo,-Co,Do);
[Acr,Bcr,Ccr,Dcr] = stabproj(Ax,Bx,Cx,Dx);
AC_1=mksys(Acr,Bcr,Ccr,Dcr,'ss');
[Acr,Bcr,Ccr,Dcr]=w_sysred(AC_1,[],[],def_flag,thr,nn);
[Aii,Bii,Cii,Dii]=ssinv(-Ai,Bi,-Ci,Di);
[Aoi,Boi,Coi,Doi]=ssinv(-Ao,Bo,-Co,Do);
[Ax,Bx,Cx,Dx]=series(Aii,Bii,Cii,Dii,Acr,Bcr,Ccr,Dcr);
[Ax,Bx,Cx,Dx]=series(Ax,Bx,Cx,Dx,Aoi,Boi,Coi,Doi);
[Acr,Bcr,Ccr,Dcr] = stabproj(Ax,Bx,Cx,Dx);
cut=length(Acr);

if length(Aca)>0
[Acr,Bcr,Ccr,Dcr] = addss(Acr,Bcr,Ccr,Dcr,Aca,Bca,Cca,Dca);
end

if cut==nc
[Acr,Bcr,Ccr,Dcr]=branch(AC_);
end
