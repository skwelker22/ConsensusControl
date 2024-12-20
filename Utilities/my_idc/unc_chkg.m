function [stchk,MUNCe,S_cl]=unc_chkg(N_P,D_P,C_i,frun,AUNC,flag,S_W);
%function [stchk,MUNCe,S_cl]=unc_chkg(N_P,D_P,C_i,frun,AUNC,flag,S_W);
% Forms the closed loop system and effective multiplicative
% uncertainty for the feedback structure.   
%                        AUNC(:,2)-->|<--AUNC(:,3)-
%                         |          |            |
%       u2---> + --> C_i --> N_P --> + --> (D_P) ---> S_W --->y1
%              |                                           |
%             -|                                           |
%              |<------------------------------------------|
%
%   N_P,D_P are packed state space descriptions of the plant
%     left factorization
%   S_W is the packed state-space representation of an output weight
%   C_i is the packed controller state space representation, possibly
%      containing an inner loop system
%   frun,AUNC(:,2,3) is the coprime factor uncertainty
% 
%  On return, stchk is the stability check, MUNCe is the 
%  bound on the effective closed loop multiplicative uncertainty
%  and S_cl is the packed state space description of u2 --> y1

%  [Adi,Bdi,Cdi,Ddi]=ssinv(A_D,B_D,-C_D,eye(size(D_D))-D_D);

if nargin<7,S_W=[];end
[Adi,Bdi,Cdi,Ddi]=branch(D_P);
[A_N,B_N,C_N,D_N]=branch(N_P);
[ninp,noutp]=size(D_N');nix=min(ninp,noutp);
[A_C,B_C,C_C,D_C]=branch(C_i);
  if ~isempty(S_W)
    [A_W,B_W,C_W,D_W]=branch(S_W);
    [A_C,B_C,C_C,D_C]=series(A_W,B_W,C_W,D_W,A_C,B_C,C_C,D_C);
    sssw=sv3_5(A_W,B_W,C_W,D_W,1,frun);rho_W=(sssw(1,:)./sssw(nix,:))';
  else
    rho_W=0*frun+1;
  end

[Apr,Bpr,Cpr,Dpr]=series(A_N,B_N,C_N,D_N,Adi,Bdi,Cdi,Ddi);
[Algn,Blgn,Clgn,Dlgn]=series(A_C,B_C,C_C,D_C,Apr,Bpr,Cpr,Dpr);
[Acly,Bcly,Ccly,Dcly]=feedbk(Algn,Blgn,Clgn,Dlgn,2);
[Acle,Bcle,Ccle,Dcle]=feedbk(Algn,Blgn,Clgn,Dlgn,1);
[Aclu,Bclu,Cclu,Dclu]=feedbk(A_C,B_C,C_C,D_C,3,Apr,Bpr,Cpr,Dpr);
[Act,Bct,Cct,Dct]=series(Acle,Bcle,Ccle,Dcle,A_C,B_C,C_C,D_C);

if size(AUNC)*[0;1]>=5
C_=ss(A_C,B_C,C_C,D_C);
P_=ss(Apr,Bpr,Cpr,Dpr);
Sy_ = fbk(P_*C_,eye(size(Dpr*D_C)));
M_ = ss(Adi,Bdi,Cdi,Ddi);
CSM_ = C_*Sy_*M_;SM_=Sy_*M_;
ssk1=sv3_5(CSM_.a,CSM_.b,CSM_.c,CSM_.d,1,frun);ssk1=ssk1(1,:)';
ssk2=sv3_5(SM_.a,SM_.b,SM_.c,SM_.d,1,frun);ssk2=ssk2(1,:)';
       loglog(frun,ssk1,frun,1./AUNC(:,4),frun,ssk2,frun,1./AUNC(:,5));grid;title('CF sensitifities and inv.uncertainty'); 
       ylabel('CSM, UNCinv, SM, UNCinv');pause
       loglog(frun,min(ssk1.*AUNC(:,4),ssk2.*AUNC(:,5)));grid;title('RSC estimate'); pause
end

% evaluate effective output multiplicative from uncertainties in sin

[Ak1,Bk1,Ck1,Dk1]=series(Adi,Bdi,Cdi,Ddi,Aclu,Bclu,Cclu,Dclu);
[Ak2,Bk2,Ck2,Dk2]=series(Adi,Bdi,Cdi,Ddi,Acle,Bcle,Ccle,Dcle);

ssk1=sv3_5(Ak1,Bk1,Ck1,Dk1,1,frun);ssk1=ssk1(1,:)';
ssk2=sv3_5(Ak2,Bk2,Ck2,Dk2,1,frun);ssk2=ssk2(1,:)';
ssk3=sv3_5(Act,Bct,Cct,Dct,1,frun);ssk3=ssk3(1,:)';
ssT=sv3_5(Acly,Bcly,Ccly,Dcly,1,frun);
%   sscs=sv3_5(Aclu,Bclu,Cclu,Dclu,1,frun);sscs=sscs(1,:)';
z_tem1=0*frun;z_tem2=z_tem1;
ind1=find(AUNC(:,2)>1.e-6);ind2=find(AUNC(:,3)>1.e-6);
z_tem1(ind1)=AUNC(ind1,2);z_tem2(ind2)=AUNC(ind2,3);
stchk=ssk1.*z_tem1+ssk2.*z_tem2;

if flag>=0
if size(AUNC)*[0;1]>=5
       loglog(frun,stchk,frun,min(ssk1.*AUNC(:,4),ssk2.*AUNC(:,5)));grid;title('Stab.Cond. check, RSC estimate'); pause
else
  loglog(frun,stchk);grid
  title('Stability Condition Check (<1)');
end
end

%------------------------------- Inner loop effective output mult.unc
MUNCe1=abs(ssk2.*(z_tem1.*ssk3./(ssT(nix,:)')+z_tem2)./(1-stchk));
MUNCe=MUNCe1.*rho_W;
if flag>=0 & nargout>1
%  loglog(frun,MUNCe,frun,MUNCe1);grid
%  title('Effective Closed-loop Multiplicative Uncertainty');
end
%----------------------------------------------------------------
if ~isempty(S_W)
  [Apr,Bpr,Cpr,Dpr]=series(A_N,B_N,C_N,D_N,Adi,Bdi,Cdi,Ddi);
  [A_C,B_C,C_C,D_C]=branch(C_i);
  [Algn,Blgn,Clgn,Dlgn]=series(A_C,B_C,C_C,D_C,Apr,Bpr,Cpr,Dpr);
  [Algn,Blgn,Clgn,Dlgn]=series(Algn,Blgn,Clgn,Dlgn,A_W,B_W,C_W,D_W);
  [Acly,Bcly,Ccly,Dcly]=feedbk(Algn,Blgn,Clgn,Dlgn,2);
end
S_cl=mksys(Acly,Bcly,Ccly,Dcly,'ss');


%disp('Nominal performance')
%  [p_a,p_b,p_c,p_d]=series(A_N,B_N,C_N,D_N,Adi,Bdi,Cdi,Ddi);
%  [c_a,c_b,c_c,c_d]=branch(C_i);
%  if ~isempty(S_W)
%     [p_a,p_b,p_c,p_d]=series(p_a,p_b,p_c,p_d,A_W,B_W,C_W,D_W);
%  end

%   scc=h_perfev(c_a,c_b,c_c,c_d,p_a,p_b,p_c,p_d,frun,2);
