function [AN,BN,CN,DN,AD,BD,CD,DD,Z]=norm_cop(AA,B,C,D,shif,lr);

%function [AN,BN,CN,DN,AD,BD,CD,DD,Z]=norm_cop(A,B,C,D,shif,lr);
% Find normalized left/right coprime factorization for a system [A,B,C,D]
% P =  DLinv NL : DL*DL' + NL*NL' = I
% P =  NR DRinv : DR'*DR + NR'*NR = I
% NL/NR = [AN,BN,CN,DN], DL/DR = [AD,BD,CD,DD]
%    shif (optional) contains the guaranteed stability margin of
%        the coprime factors; in this case, the shifted
%        DL, NL are normalized
%    lr = left = 0 / right = 1 (def=0)
% Z are the solutions to the control or filter ARE

if nargin <2; B=[];end
if nargin <5; shif=0;end
if isempty(shif);shif=0;end
if nargin <6; lr=0;end
if lr ~= 1;lr=0;end

if isempty(B), [A,B,C,D]=branch(AA);else,A=AA;end

nA=length(A);[nout,ninp]=size(D);
IA=eye(nA,nA);IO=eye(nout,nout);OI=eye(ninp,ninp);
S=OI+D'*D;Si=inv(S);Si2=sqrtm(Si);
R=IO+D*D';Ri=inv(R);Ri2=sqrtm(Ri);
AC=A-B*Si*D'*C+shif*IA;
AO=A-B*D'*Ri*C+shif*IA;

if lr ==1
   Q=C'*Ri*C;
%   [u,s,v]=svd(Q);Q=u*(s)*u';
   Z=are(AC,B*Si*B',Q);
   H=-Si*(B'*Z+D'*C);
   AN=A+B*H;CN=C+D*H;BN=B*Si2;DN=D*Si2;AD=AN;CD=H;BD=BN;DD=Si2;
else
   Q=B*Si*B';
%   [u,s,v]=svd(Q);Q=u*(s)*u';
   Z=are(AO',C'*Ri*C,Q);
   H=-(Z*C'+B*D')*Ri;
   AN=A+H*C;BN=B+H*D;CN=Ri2*C;DN=Ri2*D;AD=AN;BD=H;CD=CN;DD=Ri2;
end
if nargout <4,AN=ss(AN,BN,CN,DN);BN=ss(AD,BD,CD,DD);CN=Z;end
