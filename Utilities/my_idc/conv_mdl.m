function [A,B,C,D,F,By,Dy]=conv_mdl(At,Bt,Ct,Dt,shif);

%function [A,B,C,D,F,By,Dy,fr,MR]=mod_red(At,Bt,Ct,Dt,shif);
% model order reduction of an (A,B,C,D);
% returns the reduced order (A,B,C,D), (F,By,Dy) for its
%    nonminimal fractional representation

if nargin < 5, shif=1;end
frw=logspace(-2,2,100)';
A1=At;B1=Bt;C1=Ct;D1=Dt;
if isstr(At)
   eval(['load ',At]);A1=Apr;B1=Bpr;C1=Cpr;D1=Dpr;
end
if isstr(Bt)
   eval(['load ',Bt]);A2=Apr;B2=Bpr;C2=Cpr;D2=Dpr;
   [A1,B1,C1,D1]=series(A1,B1,C1,D1,A2,B2,C2,D2);
end
if isstr(Ct)
   eval(['load ',Ct]);A2=Apr;B2=Bpr;C2=Cpr;D2=Dpr;
   [A1,B1,C1,D1]=series(A1,B1,C1,D1,A2,B2,C2,D2);
end
if isstr(Dt)
   eval(['load ',Dt]);A2=Apr;B2=Bpr;C2=Cpr;D2=Dpr;
   [A1,B1,C1,D1]=series(A1,B1,C1,D1,A2,B2,C2,D2);
end

At=A1;Bt=B1;Ct=C1;Dt=D1;
[noutp,ninp]=size(Dt);
nA=length(At);
B=Bt;D=Dt;A=At;C=Ct;
%-------------------------- nonminimal left factorization: 
%-------------------------- Fpr,[Bpr-Bypr*Dpr,Bypr],Cpr,[Dpr,Dypr]
 while shif>=0
    IA=eye(nA,nA);IO=10*eye(noutp,noutp);
    Lpr=lqr(At'+shif*IA,Ct',IA,IO);
    F=At-Lpr'*Ct;By=Lpr';Dy=zeros(noutp,noutp);
    s1=sv3_5(F,B-By*D,C,D,1,frw);s2=sv3_5(F,-By,C,-Dy+IO/10,1,frw);
    loglog(frw,s1,'b',frw,s2,'r');title('SV of num (b) and den (r)');
    shif=input('exp.stability margin for nonminimal model= ');
    if isempty(shif);shif=-1;end 
end
Apr=A;Bpr=B;Cpr=C;Dpr=D;Fpr=F;Bypr=By;Dypr=Dy;
savsys=input(['Give filename to save ID results:   '],'s');
if length(savsys)>=2
eval(['save ',savsys,' Apr Bpr Cpr Dpr Fpr Bypr Dypr X0 DC frw sspla MR frun MUNC AUNC F_bw F_ord cl_spec P_num P_den IC_pts m_red LSMN SCALES FLAGS'])
end
disp('You NEED to re-run mimoidi to compute uncertainty estimates')
