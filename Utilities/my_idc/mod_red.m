function [A,B,C,D,F,By,Dy,fr,MR]=mod_red(At,Bt,Ct,Dt,Ft,Byt,Dyt,fil,flag,m_r);

%function [A,B,C,D,F,By,Dy,fr,MR]=mod_red(At,Bt,Ct,Dt,Ft,Byt,Dyt,fil,flag,m_r);
% model order reduction of an (A,B,C,D);
% returns the reduced order (A,B,C,D), (F,By,Dy) for its
%    nonminimal fractional representation and the frequency
%    response of the reduced and original models

A=[];B=[];C=[];D=Dt;F=[];By=[];Dy=[];fr=[];MR=[];
if nargin < 5,Ft=[];Byt=[];Dyt=[];fil=1;flag=0;end
if nargin < 8, fil=1;flag=0;end
if nargin < 9, flag=0;end
if nargin <10, m_r=[1 3 0]; end
if length(m_r)~=3;m_r=[1 3 0];end
if flag > 2; m_r(2)=3;end
if m_r(2)==1 & m_r(3)<0;m_r(3)=1.e-4;end
if m_r(2)==2 & m_r(3)>0;m_r(3)=-40;end

[noutp,ninp]=size(Dt);nA=length(At);


if m_r(1)==-1
   [At,Bt,Ct,Dt]=ssinv(At,Bt,Ct,Dt);
end

    [A,B,C,D]=sys_red(At,Bt,Ct,Dt,3,[0 1 0]);
     nA=length(A);
%-------------------------- nonminimal left factorization: 
%-------------------------- Fpr,[Bpr-Bypr*Dpr,Bypr],Cpr,[Dpr,Dypr]
     shif=mean(fil)/10;
%     shif=input('exp.stability margin for nonminimal model= ');
%     if length(shif)<1; shif=mean(fil)/2,end
     IA=eye(nA,nA);IO=10*eye(noutp,noutp);
     Lpr=lqr(A'+shif*IA,C',IA,IO);
     F=A-Lpr'*C;By=Lpr';Dy=zeros(noutp,noutp);

