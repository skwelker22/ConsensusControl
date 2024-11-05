function [stchk,frun,MUNCe,S_cl]=unc_chk(s_in,A,B,C,D,u,y,t,flag);
%function [stchk,frun,MUNCe,S_cl]=unc_chk(s_in,A,B,C,D,u,y,t,flag);
% Forms the inner loop system and outer loop effective plant for the
% feedback structure
%       u2---> + ---> C1_ ---> S1_ ------>y1
%              |                    |
%             -|                    |
%              |<-------------------|
% s_in is an id file containing the S1_ model (Apr,Bpr,...,Fpr,...) and 
%      uncertainty info (frun,MUNC,AUNC)
% A,B,C,D is the controller state-space realization
% u,y,t is the optional i/o info for re-calclulation of 
%      coprime factor uncertainty decomposition
% plots stchk vs. frun if flag>=0

if nargin < 9; flag=0;end
if nargin < 6; u=[];end
flagt=flag;
eval(['load ' s_in]);
[ninp,noutp]=size(Dpr');nix=min(ninp,noutp);

[Algn,Blgn,Clgn,Dlgn]=series(A,B,C,D,Apr,Bpr,Cpr,Dpr);
[Acly,Bcly,Ccly,Dcly]=feedbk(Algn,Blgn,Clgn,Dlgn,2);
[Acle,Bcle,Ccle,Dcle]=feedbk(Algn,Blgn,Clgn,Dlgn,1);
[Aclu,Bclu,Cclu,Dclu]=feedbk(A,B,C,D,3,Apr,Bpr,Cpr,Dpr);
[Act,Bct,Cct,Dct]=series(Acle,Bcle,Ccle,Dcle,A,B,C,D);

% evaluate effective output multiplicative from uncertainties in sin

[Adi,Bdi,Cdi,Ddi]=ssinv(Fpr,Bypr,-Cpr,eye(size(Dypr))-Dypr);
[Ak1,Bk1,Ck1,Dk1]=series(Adi,Bdi,Cdi,Ddi,Aclu,Bclu,Cclu,Dclu);
[Ak2,Bk2,Ck2,Dk2]=series(Adi,Bdi,Cdi,Ddi,Acle,Bcle,Ccle,Dcle);
ssk1=sv3_5(Ak1,Bk1,Ck1,Dk1,1,frun);ssk1=ssk1(1,:)';
ssk2=sv3_5(Ak2,Bk2,Ck2,Dk2,1,frun);ssk2=ssk2(1,:)';
ssk3=sv3_5(Act,Bct,Cct,Dct,1,frun);ssk3=ssk3(1,:)';

if length(u)<2
stchk=ssk1.*AUNC(:,2)+ssk2.*AUNC(:,3);
else
% ssT=sv3_5(Acly,Bcly,Ccly,Dcly,1,frun);ssT=ssT(1,:)';
% ssS=sv3_5(Acle,Bcle,Ccle,Dcle,1,frun);ssS=ssS(1,:)';
stp=t(2)-t(1);

   nan_ind=[find(isnan(u(:,1)));nn+1];
   n_nan=length(nan_ind);i_nnan=find(~isnan(u(:,1)));

X0=[];     i_nan1=1;
for i_nan=1:length(nan_ind);
   ind_nan=[i_nan1:nan_ind(i_nan)-1];
   i_nan1=nan_ind(i_nan)+1;
   UY_nan=[u(ind_nan,:),y(ind_nan,:)];
   Y_nan=y(ind_nan,:);T_nan=t(ind_nan);T_nan=T_nan-T_nan(1);
   xx0=icestim(Fpr,[Bpr-Bypr*Dpr,Bypr],Cpr,[Dpr,Dypr],UY_nan,Y_nan,T_nan);
   X0=[X0,xx0];
end

wu=lsim_nan(Fpr,Bpr-Bypr*Dpr,Cpr,Dpr,u,t,X0);
wy=lsim_nan(Fpr,Bypr,Cpr,Dypr,y,t);
Ferr=y-(wu+wy);
ssva=[ssk1';ssk2'];
[frun,s_e,s_u,s_y]=unc_fft(Ferr,u,y,frun,stp);
UNCFo=unc_bn(s_e,s_u,s_y,ssva,2);
% stchk=unc_chk0(ssS,ssT,Apr,Fpr,Bpr,Cpr,Dpr,Bypr,Dypr,u,y,t,frun,stp,flag);
AUNC(:,2:3)=UNCFo;
z_tem1=0*frun;z_tem2=z_tem1;
ind1=find(AUNC(:,2)>1.e-6);ind2=find(AUNC(:,3)>1.e-6);
z_tem1(ind1)=AUNC(ind1,2);z_tem2(ind2)=AUNC(ind2,3);
stchk=ssk1.*z_tem1+ssk2.*z_tem2;
end

if flag>=0
  disp('Nominal Inner loop eigenvalues'),
     if flag>2;eig(Acle);else;max(real(eig(Acle))),end
  loglog(frun,stchk);grid
  title('Stability Condition Check (<1)');
  pause(1)
end

ssT=sv3_5(Acly,Bcly,Ccly,Dcly,1,frun);
MUNCe=ssk2.*(z_tem1.*ssk3./(ssT(nix,:)')+z_tem2)./(1-stchk);
S_cl=mksys(Acly,Bcly,Ccly,Dcly,'ss');
