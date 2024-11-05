function [AW,BW,CW,DW]=sp_facl(A1,B1,C1,D1,A2,B2,C2,D2,tol);

%function [AW,BW,CW,DW]=sp_facl(A1,B1,C1,D1,A2,B2,C2,D2,tol);
%   Given a factorization G1 G2inv find a spectral factor W
%   such that W~W = G~G where G = [G1 ; G2].
%   For G2 = I, give A2=[].

if nargin < 8; A2=[]; end
if nargin < 9; tol=1e-6;  end

if ~isempty(A2)
   [A,B,C,D] = append(A1,B1,C1,D1,A2,B2,C2,D2)
   I2=eye(size(D2));
   B=B*[I2;I2];  D=D*[I2;I2];
else
   A=A1;B=B1;C=C1;D=D1;
end

[UD,SD,VD]=svd(D);  [rd,cd]=size(D); nsd1=min(rd,cd);
SD1=diag(max(diag(SD),tol));
SD(1:nsd1,1:nsd1)=SD1(1:nsd1,1:nsd1);
D=UD*SD*VD';

R=D'*D;
Ri=inv(R);Ri2=sqrtm(Ri);R2=sqrtm(R);

AY=(A-B*Ri*D'*C);
BY=B*Ri*B';
I1=eye(size(D*D'));
CY=C'*(I1-D*Ri*D')*C;
Y=are(AY,BY,CY);

AW=A; BW=B; DW=R2; CW=Ri2*(D'*C+B'*Y);

