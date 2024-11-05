function [D_N,D_D,Q_]=ncf_dif2(P_2,P_1N,P_1D,shift);

%function [D_N,D_D,Q_]=ncf_dif2(P_2,P_1N,P_1D,shift);
%   Given two systems, in packed format, find an optimal NCF error
%   1 = nominal system, 2 = perturbation.
%   P_1 can be defined in terms of a coprime factorization N and D
%     otherwise an NCF is computed.
%   On return, D_N and D_D are the optimal uncertainty factors
%     and Q_ is the optimal weight (recovered by branch)

if nargin<4;shift=0;end

if nargin<3 | isempty(P_1D)
   [A1,B1,C1,D1]=branch(P_1N);
   [ANL1,BNL1,CNL1,DNL1,ADL1,BDL1,CDL1,DDL1]=norm_cop(A1,B1,C1,D1,shift,0);
else
   [ANL1,BNL1,CNL1,DNL1]=branch(P_1N);
   [ADL1,BDL1,CDL1,DDL1]=branch(P_1D);
end

[A2,B2,C2,D2]=branch(P_2);
[ANL2,BNL2,CNL2,DNL2,ADL2,BDL2,CDL2,DDL2]=norm_cop(A2,B2,C2,D2,shift,0);
[ANR2,BNR2,CNR2,DNR2,ADR2,BDR2,CDR2,DDR2]=norm_cop(A2,B2,C2,D2,shift,1);

%[AC1,BC1,CC1,DC1]=series(ADR2,BDR2,CDR2,DDR2,ANL1,BNL1,CNL1,DNL1);
%[AC2,BC2,CC2,DC2]=series(ANR2,BNR2,CNR2,DNR2,ADL1,BDL1,CDL1,DDL1);
%[AC,BC,CC,DC]=addss(AC1,-BC1,CC1,-DC1,AC2,BC2,CC2,DC2);

[AC1,BC1,CC1,DC1]=series(-ANL2',-CNL2',BNL2',DNL2',ANL1,BNL1,CNL1,DNL1);
[AC2,BC2,CC2,DC2]=series(-ADL2',-CDL2',BDL2',DDL2',ADL1,BDL1,CDL1,DDL1);
[AA,BA,CA,DA]=addss(AC1,BC1,CC1,DC1,AC2,BC2,CC2,DC2);

disp('ncf-diff, Solving Nehari...')
%[AQ,BQ,CQ,DQ,sigr]=nehari(AA,BA,CA,DA);
%disp(['ncf-diff, Nehari sigma = ',num2str(sigr)])
[AQ,BQ,CQ,DQ,sigr]=stabproj(AA,BA,CA,DA);DQ=DA;
[AC1,BC1,CC1,DC1]=series(ANL2,BNL2,CNL2,DNL2,AQ,BQ,CQ,DQ);
[AC2,BC2,CC2,DC2]=series(ADL2,BDL2,CDL2,DDL2,AQ,BQ,CQ,DQ);

[ADN,BDN,CDN,DDN]=addss(ANL1,BNL1,CNL1,DNL1,AC1,-BC1,CC1,-DC1);
[ADD,BDD,CDD,DDD]=addss(ADL1,BDL1,CDL1,DDL1,AC2,-BC2,CC2,-DC2);

%Q_=mksys(AQ,BQ,CQ,DQ,'ss');
%D_N=mksys(ADN,BDN,CDN,DDN,'ss');
%D_D=mksys(ADD,BDD,CDD,DDD,'ss');

Q_=ss(AQ,BQ,CQ,DQ);
D_N=ss(ADN,BDN,CDN,DDN);
D_D=ss(ADD,BDD,CDD,DDD);
