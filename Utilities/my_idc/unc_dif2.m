function [D_N,D_D,ssn,ssd,sc,ssw]=unc_diff(P_,P0_N,P0_D,C_,W_,fr,tols);

%function [D_N,D_D,ssn,ssd,sc,Q_]=unc_diff(P_,P0_N,P0_D,C_,W_,fr,tols);
% Compute CF error bounds.  
% P_ is the packed state-space representation of the perturbed system
% P0_N, P0_D are the packed coprime factors of the nominal system; if P0_D=[]
%   then an NCF of P0_N is used. Alternatively, P0_N can be a string specifying 
%   the name of an id-file for P0.
% C_ is the packed controller state-space or the name of a controller-file
% fr = frequency vector for sv plots
% tols = tolerance vector [D_level,e_shift,red_tol]


if nargin<5;W_=[];end
if nargin<6;fr=logspace(-3,3,100)';end
if isempty(fr);fr=logspace(-3,3,100)';end
if nargin<7;tols=[1e-5,1e-4,1e-5];end
if length(tols)<2;tols(2)=1e-4;end
if length(tols)<3;tols(3)=1e-5;end
red_tol=tols(3);shift=tols(2);d_lev=tols(1);

[Ap,Bp,Cp,Dp]=branch(P_);Ap=sh_m(Ap,shift);Dp=fix_d(Dp,d_lev);
if isstr(P0_N)
   eval(['load ',P0_N])
   Dpr=fix_d(Dpr,d_lev);Apr=sh_m(Apr,shift);
   AD0=sh_m(Fpr,shift);BD0=Bypr;CD0=-Cpr;DD0=eye(size(Dypr))-Dypr;
   AN0=sh_m(Fpr,shift);BN0=Bpr-Bypr*Dpr;CN0=Cpr;DN0=Dpr;
   [AN,BN,CN,DN,AD,BD,CD,DD]=norm_cop(Apr,Bpr,Cpr,Dpr,0,0);
   P_2=mksys(Apr,Bpr,Cpr,Dpr,'ss');
   P_1N=mksys(AN0,BN0,CN0,DN0,'ss');
   P_1D=mksys(AD0,BD0,CD0,DD0,'ss');
   [D_Nb,D_Db,Q_0]=ncf_dif2(P_2,P_1N,P_1D,0);
   [Aq,Bq,Cq,Dq]=branch(Q_0);
   ssq=sv3_5(Aq,Bq,Cq,Dq,1,fr);
   TZ=tzero(Aq,Bq,Cq,Dq);max(real(TZ))
   clf
   loglog(fr,ssq);title('Q');pause
else
   if isempty(P0_D)
      [a,b,c,d]=branch(P0_N);a=sh_m(a,shift);d=fix_d(d,d_lev);
      [AN,BN,CN,DN,AD,BD,CD,DD]=norm_cop(a,b,c,d,0,0);
   else
      [AN,BN,CN,DN]=branch(P0_N);AN=sh_m(AN,shift);DN=fix_d(DN,d_lev);
      [AD,BD,CD,DD]=branch(P0_D);AD=sh_m(AD,shift);
   end
end
[AN,BN,CN]=obalreal(AN,BN,CN);
[AD,BD,CD]=obalreal(AD,BD,CD);

if isstr(C_)
   eval(['load ',C_])
else
   [Acr,Bcr,Ccr,Dcr]=branch(C_);
end
Acr=sh_m(Acr,shift);
Dcr=fix_d(Dcr,d_lev);

[ADi,BDi,CDi,DDi]=ssinv(AD,BD,CD,DD);
[ADi0,BDi0,CDi0,DDi0]=ssinv(AD0,BD0,CD0,DD0);

%--- e->y system
[Aey,Bey,Cey,Dey]=series(Acr,Bcr,Ccr,Dcr,AN,BN,CN,DN);
[Aey,Bey,Cey,Dey]=feedbk(ADi,BDi,CDi,DDi,3,Aey,Bey,Cey,Dey);
[Aey0,Bey0,Cey0,Dey0]=series(Acr,Bcr,Ccr,Dcr,AN0,BN0,CN0,DN0);
[Aey0,Bey0,Cey0,Dey0]=feedbk(ADi0,BDi0,CDi0,DDi0,3,Aey0,Bey0,Cey0,Dey0);
%--- e->u system
[Aeu,Beu,Ceu,Deu]=series(ADi,BDi,CDi,DDi,Acr,-Bcr,Ccr,-Dcr);
[Aeu,Beu,Ceu,Deu]=feedbk(Aeu,Beu,Ceu,Deu,3,AN,-BN,CN,-DN);
[Aeu0,Beu0,Ceu0,Deu0]=series(ADi0,BDi0,CDi0,DDi0,Acr,-Bcr,Ccr,-Dcr);
[Aeu0,Beu0,Ceu0,Deu0]=feedbk(Aeu0,Beu0,Ceu0,Deu0,3,AN0,-BN0,CN0,-DN0);

[Aey,Bey,Cey]=obalreal(Aey,Bey,Cey);
[Aeu,Beu,Ceu]=obalreal(Aeu,Beu,Ceu);
[Aey0,Bey0,Cey0]=obalreal(Aey0,Bey0,Cey0);
[Aeu0,Beu0,Ceu0]=obalreal(Aeu0,Beu0,Ceu0);

AY_=mksys(Aey,Bey,Cey,Dey,'ss');
AU_=mksys(Aeu,Beu,Ceu,Deu,'ss');
[Aey,Bey,Cey,Dey]=w_sysred(AY_,[],[],0,0,[3,red_tol]);
[Aeu,Beu,Ceu,Deu]=w_sysred(AU_,[],[],0,0,[3,red_tol]);

%--- spectral factors
[AAi,BAi,CAi,DAi]=sp_facr(Aeu,Beu,Ceu,Deu);
[ABi,BBi,CBi,DBi]=sp_facr(Aey,Bey,Cey,Dey);

done=1;
while done~=0
   
   if ~isempty(W_)
      BAi=BAi*W_;DAi=DAi*W_;  
   end
   STEM=mksys(AAi,BAi,CAi,DAi,'ss');[AAi,BAi,CAi,DAi]=w_sysred(STEM,[],[],0,0,[3,red_tol]);
   [AA,BA,CA,DA]=ssinv(AAi,BAi,CAi,DAi);
   [AB,BB,CB,DB]=ssinv(ABi,BBi,CBi,DBi);
      
   [ANb,BNb,CNb,DNb]=series(AAi,BAi,CAi,DAi,AN,BN,CN,DN);
   [ADb,BDb,CDb,DDb]=series(ABi,BBi,CBi,DBi,AD,BD,CD,DD);
   [Apb,Bpb,Cpb,Dpb]=series(AAi,BAi,CAi,DAi,Ap,Bp,Cp,Dp);
   [Apb,Bpb,Cpb,Dpb]=series(Apb,Bpb,Cpb,Dpb,AB,BB,CB,DB);
   STEM=mksys(ANb,BNb,CNb,DNb,'ss');[ANb,BNb,CNb,DNb]=w_sysred(STEM,[],[],0,0,[3,red_tol]);
   STEM=mksys(ADb,BDb,CDb,DDb,'ss');[ADb,BDb,CDb,DDb]=w_sysred(STEM,[],[],0,0,[3,red_tol]);
   STEM=mksys(Apb,Bpb,Cpb,Dpb,'ss');[Apb,Bpb,Cpb,Dpb]=w_sysred(STEM,[],[],0,0,[3,red_tol]);
   
   P_2=mksys(Apb,Bpb,Cpb,Dpb,'ss');
   P_1N=mksys(ANb,BNb,CNb,DNb,'ss');
   P_1D=mksys(ADb,BDb,CDb,DDb,'ss');
   
   [D_Nb,D_Db,Q_]=ncf_dif2(P_2,P_1N,P_1D,0);
   
   [ADNb,BDNb,CDNb,DDNb]=branch(D_Nb);
   [ADDb,BDDb,CDDb,DDDb]=branch(D_Db);
   [ADN,BDN,CDN,DDN]=series(AA,BA,CA,DA,ADNb,BDNb,CDNb,DDNb);
   [ADD,BDD,CDD,DDD]=series(AB,BB,CB,DB,ADDb,BDDb,CDDb,DDDb);
   [Aqi,Bqi,Cqi,Dqi]=ssinv(Aq,Bq,Cq,Dq);
   [ADN,BDN,CDN,DDN]=series(ADN,BDN,CDN,DDN,Aqi,Bqi,Cqi,Dqi);
   [ADD,BDD,CDD,DDD]=series(ADD,BDD,CDD,DDD,Aqi,Bqi,Cqi,Dqi);
   
   %   [ADN,BDN,CDN,DDN]=branch(D_Nb);
   %   [ADD,BDD,CDD,DDD]=branch(D_Db);
   
   ADN=sh_m(ADN,-shift);ADD=sh_m(ADD,-shift);
   D_N=mksys(ADN,BDN,CDN,DDN,'ss');
   D_D=mksys(ADD,BDD,CDD,DDD,'ss');
%   [ADN,BDN,CDN,DDN]=w_sysred(D_N,[],[],0,0,[3,red_tol]);
%   [ADD,BDD,CDD,DDD]=w_sysred(D_D,[],[],0,0,[3,red_tol]);
   D_N=mksys(ADN,BDN,CDN,DDN,'ss');
   D_D=mksys(ADD,BDD,CDD,DDD,'ss');
   
   sscsd=sv3_5(sh_m(Aeu0,-shift),Beu0,Ceu0,Deu0,1,fr);
   sssd=sv3_5(sh_m(Aey0,-shift),Bey0,Cey0,Dey0,1,fr);
   ssn=sv3_5(ADN,BDN,CDN,DDN,1,fr);
   ssd=sv3_5(ADD,BDD,CDD,DDD,1,fr);
%   ssnb=sv3_5(ADNb,BDNb,CDNb,DDNb,1,fr);
%   ssdb=sv3_5(ADDb,BDDb,CDDb,DDDb,1,fr);
   ssw=sqrt(ssd(1,:)./ssn(1,:).*sscsd(1,:)./sssd(1,:));
   sc1=sscsd(1,:).*ssn(1,:);sc2=sssd(1,:).*ssd(1,:);
   sc=sc1+sc2;
   
   loglog(fr,sc,fr,sc1,fr,sc2);title('RSC');pause
   loglog(fr,sqrt(ssw));title('w-opt');pause
   loglog(fr,ssn,'r',fr,ssd,'b');title('Delta_N (r) and Delta_D (b) SV');pause
   %loglog(fr,ssnb,'r',fr,ssdb,'b');title('Delta_Nb (r) and Delta_Db (b) SV');pause
   
   %   W_=u_d_app(fr,ssw',length(DAi));
   
   done=input('done? [0=y]  ');
   if isempty(done);done=0;end
   W_=done;
end

