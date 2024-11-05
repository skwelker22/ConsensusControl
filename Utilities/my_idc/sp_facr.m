function [AW,BW,CW,DW]=sp_facr(A1,B1,C1,D1,A2,B2,C2,D2,tol);

%function [AW,BW,CW,DW]=sp_facr(A1,B1,C1,D1,A2,B2,C2,D2,tol);
%   Given a factorization G2inv G1 find a spectral factor W
%   such that WW~ = GG~ where G = [G1 G2].
%   For G2 = I, give A2=[].

if nargin < 8; A2=[]; end
if nargin < 9; tol=1e-6;  end

if ~isempty(A2)
   [A,B,C,D] = append(A1,B1,C1,D1,A2,B2,C2,D2)
   I2=eye(size(D2));
   C=[I2,I2]*C;  D=[I2,I2]*D;
else
   A=A1;B=B1;C=C1;D=D1;
end

D=fix_d(D,tol);

R=D*D';
Ri=inv(R);Ri2=sqrtm(Ri);R2=sqrtm(R);

AY=(A-B*D'*Ri*C)';
BY=C'*Ri*C;
I1=eye(size(D'*D));
CY=B*(I1-D'*Ri*D)*B';
Y=are(AY,BY,CY);

AW=A; CW=C; DW=R2; BW=(B*D'+Y*C')*Ri2;

