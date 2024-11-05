function [AN,BN,CN,DN,AD,BD,CD,DD,Z]=d_ncf(A,B,C,D,shif);

%function [AN,BN,CN,DN,AD,BD,CD,DD,Z]=d_ncf(A,B,C,D,shif);
% Find normalized left coprime factorization for a discrete system [A,B,C,D]
% P =  DLinv NL : DL*DL'+NL*NL'=I
% NL=[AN,BN,CN,DN], DL=[AD,BD,CD,DD]
%    shif (optional) contains the guaranteed stability margin of
%        the coprime factors; in this case, the shifted
%        DL, NL are normalized
% Z are the solutions to the control or filter ARE

if nargin <5; shif=0;end

nA=length(A);[nout,ninp]=size(D);
IA=eye(nA,nA);IO=eye(nout,nout);OI=eye(ninp,ninp);
% S=OI+D'*D;Si=inv(S);Si2=sqrtm(Si);
% AC=A-B*Si*D'*C+shif*IA;
R=IO+D*D';Ri=inv(R);Ri2=sqrtm(Ri);
AO=(A-B*D'*Ri*C)*(1+shif);
[Lpr,Z]=dlqr(AO',C',B*(OI-D'*Ri*D)*B',Ri);
% [Kpr,Z]=lqr(AC,B,C'*(IO-D*Si*D')*C,Si);
H=-(Z*C'+B*D')*Ri;
AN=A+H*C;BN=B+H*D;CN=Ri2*C;DN=Ri2*D;AD=A+H*C;BD=H;CD=CN;DD=Ri2;
