function [N,M]=ncf(sys,lr,pc,shift);

%function [N,M]=ncf(sys,lr,pc,shif);
% Find normalized left/right coprime factorization for a ss system [A,B,C,D]
% P =  DLinv NL : DL*DL' + NL*NL' = I
% P =  NR DRinv : DR'*DR + NR'*NR = I
%    lr = left = 0 / right = 1 (def=0)
%    pc = plant 0 / controller = 1 type
%    shif (optional) contains the guaranteed stability margin of
%        the coprime factors; in this case, the shifted
%        DL, NL are normalized

if nargin <4; shift=0;end
if nargin <3; pc=0;end
if nargin <2; lr=0;end

if isobject(sys)
    [A,B,C,D]=ssdata(sys);
else
    [A,B,C,D]=branch(sys);
end
nA=length(A);[nout,ninp]=size(D);
IA=eye(nA,nA);IO=eye(nout,nout);OI=eye(ninp,ninp);
S=OI+D'*D;Si=inv(S);Si2=sqrtm(Si);
R=IO+D*D';Ri=inv(R);Ri2=sqrtm(Ri);
AC=A-B*Si*D'*C+shift*IA;
AO=A-B*D'*Ri*C+shift*IA;

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
N=ss(AN,BN,CN,DN);M=ss(AD,BD,CD,DD);
if lr ==1, G=[N;M]; K=[M;N]; else, G=[-M N]; K=[-N M];end

if nargout <2, if pc ==0, N=G; else, N=K; end, end 
