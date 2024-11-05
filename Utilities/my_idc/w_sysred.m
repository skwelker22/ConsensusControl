function [Acr,Bcr,Ccr,Dcr]=w_sysred(AC_,AWi_,AWo_,def_flag,thr,nn);
%function [Acr,Bcr,Ccr,Dcr]=w_sysred(AC_,AWi_,AWo_,def_flag,thr,nn);
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
if isempty(AWi_)
    iflag=1;Ai=[];Bi=[];Ci=[];Di=eye(size(Dc'*Dc));
else
    iflag=0;[Ai,Bi,Ci,Di]=branch(AWi_);
end
if isempty(AWo_)
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
if isempty(Ai)
    Ap=Ac;Bp=Bc;
else
    Ap=[Ac,Bc*Ci;0*Bi*Cc,Ai];Bp=[Bc*Di;Bi];
end
if isempty(Ao)
    Aq=Ac;Cq=Cc';
else
    Aq=[Ac,0*Bc*Co;Bo*Cc,Ao];Cq=[Cc'*Do';Co'];
end

PP=lyap(Ap,Bp*Bp');QQ=lyap(Aq',Cq*Cq');
P=PP(1:nc,1:nc);Q=QQ(1:nc,1:nc);

[U,S,V]=svd(P);Sp=sqrt(S);R=Sp*U';
[V,S,X]=svd(R*Q*R');Sq=sqrt(S);T=sqrt(Sq)*V'*pinv(Sp)*U';
cut=[];
dsq=diag(Sq);

if nn(1)==1
    cut=nn(2);
    disp(['Weighted Reduction cut = ',num2str(cut)])
elseif nn(1)==2
    cut=max(find(dsq>nn(2)));
    if isempty(cut);cut=0;end
    disp(['Weighted Reduction cut = ',num2str(cut)])
elseif nn(1)==3
    cut=max(find(dsq>nn(2)*max(dsq)));
    disp(['Weighted Reduction cut = ',num2str(cut)])
else
    disp('Weighted Hankel Singular values')
    if def_flag >=2;dsq',end
    semilogy(diag(Sq),'x');grid;title('Weighted Hankel SV')
    cut=input('give number of states in the reduced system [all]  ');
end

if length(cut)<1;cut=length(Sq);end
cut=max(min(cut,nc),0);

if cut>0 & nn(1) ~= 0;
    ddsq=[dsq;-1];
    while ddsq(cut)< 1.2*ddsq(cut+1) & cut<nc;
        cut=cut+1;
    end
end
if def_flag>0
    disp(['Adjusted Weighted Reduction cut = ',num2str(cut)])
end

Ti=pinv(T);
A1=T*Ac*Ti;B1=T*Bc;C1=Cc*Ti;
Acr=A1(1:cut,1:cut);Bcr=B1(1:cut,:);Ccr=C1(:,1:cut);Dcr=Dc0;

if length(Aca)>0
    [Acr,Bcr,Ccr,Dcr] = addss(Acr,Bcr,Ccr,Dcr,Aca,Bca,Cca,Dca);
end

if cut==nc
    [Acr,Bcr,Ccr,Dcr]=branch(AC_);
end

if nargout<4, Acr = ss(Acr,Bcr,Ccr,Dcr);end
