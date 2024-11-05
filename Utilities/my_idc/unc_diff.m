function [D_N,D_D,D_ND,D_DN,SVS]=unc_diff(P_,P0_N,P0_D,C_,fr,tols);

%function [D_N,D_D,D_ND,D_DN,SVS]=unc_diff(P_,P0_N,P0_D,C_,fr,tols);
% Compute CF error bounds.  
% P_ is the packed state-space representation of the perturbed system
% P0_N, P0_D are the packed coprime factors of the nominal system; if P0_D=[]
%   then an NCF of P0_N is used. Alternatively, P0_N can be a string specifying 
%   the name of an id-file for P0.
% C_ is the packed controller state-space or the name of a controller-file
% fr = frequency vector for sv plots
% tols = tolerance vector [D_level,e_shift,red_tol]

hold off;clf
if nargin<5;fr=logspace(-3,3,100)';end
if isempty(fr);fr=logspace(-3,3,100)';end
if nargin<6;tols=[1e-5,1e-4,1e-5];end
if length(tols)<2;tols(2)=1e-4;end
if length(tols)<3;tols(3)=1e-5;end
red_tol=tols(3);shift=tols(2);d_lev=tols(1);

[Ap,Bp,Cp,Dp]=branch(P_);Ap=sh_m(Ap,shift);Dp=fix_d(Dp,d_lev);
[ANX,BNX,CNX,DNX,ADX,BDX,CDX,DDX]=norm_cop(Ap,Bp,Cp,Dp,0,0);
if isstr(P0_N)
   eval(['load ',P0_N])
   Dpr=fix_d(Dpr,d_lev);Apr=sh_m(Apr,shift);
   AD0=sh_m(Fpr,shift);BD0=Bypr;CD0=-Cpr;DD0=eye(size(Dypr))-Dypr;
   AN0=sh_m(Fpr,shift);BN0=Bpr-Bypr*Dpr;CN0=Cpr;DN0=Dpr;
else
   if isempty(P0_D)
      [a,b,c,d]=branch(P0_N);a=sh_m(a,shift);d=fix_d(d,d_lev);
   else
      [AN0,BN0,CN0,DN0]=branch(P0_N);
      AN0=sh_m(AN0,shift);DN0=fix_d(DN0,d_lev);
      [AD0,BD0,CD0,DD0]=branch(P0_D);AD0=sh_m(AD0,shift);
   end
end
[AN0,BN0,CN0]=obalreal(AN0,BN0,CN0);
[AD0,BD0,CD0]=obalreal(AD0,BD0,CD0);

if isstr(C_)
   eval(['load ',C_])
else
   [Acr,Bcr,Ccr,Dcr]=branch(C_);
end
Acr=sh_m(Acr,shift);Dcr=fix_d(Dcr,d_lev);

[ADi0,BDi0,CDi0,DDi0]=ssinv(AD0,BD0,CD0,DD0);

%--- e->y system
[Aey0,Bey0,Cey0,Dey0]=series(Acr,Bcr,Ccr,Dcr,AN0,BN0,CN0,DN0);
[Aey0,Bey0,Cey0,Dey0]=feedbk(ADi0,BDi0,CDi0,DDi0,3,Aey0,Bey0,Cey0,Dey0);
%--- e->u system
[Aeu0,Beu0,Ceu0,Deu0]=series(ADi0,BDi0,CDi0,DDi0,Acr,-Bcr,Ccr,-Dcr);
[Aeu0,Beu0,Ceu0,Deu0]=feedbk(Aeu0,Beu0,Ceu0,Deu0,3,AN0,-BN0,CN0,-DN0);
%[Aey0,Bey0,Cey0]=obalreal(Aey0,Bey0,Cey0);
%[Aeu0,Beu0,Ceu0]=obalreal(Aeu0,Beu0,Ceu0);

AY_=mksys(Aey0,Bey0,Cey0,Dey0,'ss');
AU_=mksys(Aeu0,Beu0,Ceu0,Deu0,'ss');
[Aey0,Bey0,Cey0,Dey0]=w_sysred(AY_,[],[],0,0,[3,red_tol]);
[Aeu0,Beu0,Ceu0,Deu0]=w_sysred(AU_,[],[],0,0,[3,red_tol]);

%--- spectral factors
[AAi,BAi,CAi,DAi]=sp_facr(Aeu0,Beu0,Ceu0,Deu0);
[ABi,BBi,CBi,DBi]=sp_facr(Aey0,Bey0,Cey0,Dey0);
[AA,BA,CA,DA]=ssinv(AAi,BAi,CAi,DAi);
[AB,BB,CB,DB]=ssinv(ABi,BBi,CBi,DBi);

SW_N=mksys(AAi,BAi,CAi,DAi,'ss');
SW_D=mksys(ABi,BBi,CBi,DBi,'ss');
SW_Ni=mksys(AA,BA,CA,DA,'ss');
SW_Di=mksys(AB,BB,CB,DB,'ss');
N_2=mksys(AN0,BN0,CN0,DN0,'ss');
D_2=mksys(AD0,BD0,CD0,DD0,'ss');
N_1=mksys(ANX,BNX,CNX,DNX,'ss');
D_1=mksys(ADX,BDX,CDX,DDX,'ss');

[D_S,ND_S]=ncf_diff(N_1,N_2,D_1,D_2,[],red_tol);
[N_S,DN_S]=ncf_diff(D_1,D_2,N_1,N_2,[],red_tol);

[ADN,BDN,CDN,DDN]=branch(N_S);
[ADD,BDD,CDD,DDD]=branch(D_S);
[ADND,BDND,CDND,DDND]=branch(ND_S);
[ADDN,BDDN,CDDN,DDDN]=branch(DN_S);
ADN=sh_m(ADN,-shift);ADD=sh_m(ADD,-shift);
ADND=sh_m(ADND,-shift);ADDN=sh_m(ADDN,-shift);
D_N=mksys(ADN,BDN,CDN,DDN,'ss');
D_D=mksys(ADD,BDD,CDD,DDD,'ss');
D_ND=mksys(ADND,BDND,CDND,DDND,'ss');
D_DN=mksys(ADDN,BDDN,CDDN,DDDN,'ss');

sscsd=sv3_5(sh_m(Aeu0,-shift),Beu0,Ceu0,Deu0,1,fr);
sssd=sv3_5(sh_m(Aey0,-shift),Bey0,Cey0,Dey0,1,fr);
ssn=sv3_5(ADN,BDN,CDN,DDN,1,fr);
ssd=sv3_5(ADD,BDD,CDD,DDD,1,fr);
ssnd=sv3_5(ADND,BDND,CDND,DDND,1,fr);
ssdn=sv3_5(ADDN,BDDN,CDDN,DDDN,1,fr);
ssw=sqrt(ssd(1,:)./ssn(1,:).*sscsd(1,:)./sssd(1,:));
sc1n=sscsd(1,:).*ssn(1,:);
sc1d=sssd(1,:).*ssdn(1,:);
sc2d=sssd(1,:).*ssd(1,:);
sc2n=sscsd(1,:).*ssnd(1,:);
sc1=sc1n+sc1d;sc2=sc2n+sc2d;sc=min(sc1,sc2);
disp(['max RSC = ',num2str(max(sc))]);
loglog(fr,sat([sc1n;sc1d;sc2n;sc2d],1e5,1e-5),fr,sc,'x');
title('RSC');grid;pause
loglog(fr,sqrt(ssw));title('w-opt');pause
loglog(fr,ssn,'r',fr,ssd,'b');title('Delta_N (r) and Delta_D (b) SV');pause
%loglog(fr,ssnb,'r',fr,ssdb,'b');title('Delta_Nb (r) and Delta_Db (b) SV');pause
SVS=[sc;sc1;sc2;ssn;ssd;ssnd;ssdn];