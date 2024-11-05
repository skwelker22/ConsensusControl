function [AX,BX,CX,DX,sigr]=nehari(AA,B,C,D,tol,Y);
% [AX,BX,CX,DX,sigr]=nehari(A,B,C,D,tol);
% Solves the Nehari approximation problem
%      min_X ||R - X||_inf
%        s.t., R in L_inf, X in H_inf 
%   R is in [A,B,C,D]
%   X is in [AX,BX,CX,DX]; sigr is the minimum distance
%   tol = Hankel inversion tolerance; high frequency modes are
%         cut beyond 0.1/tol. To bypass this reduction, give tol as
%         a negative number

nosf=0;
if nargin < 5; tol=1.e-8;        end
if tol==0;     tol=1e-8;         end
if tol<0;      tol=-tol; nosf=1; end
if nargin < 6; Y=0;              end

temp=version;
if str2num(temp(1))>=6;
    if isobject(AA);AA=ss(AA);A=AA.a;B=AA.b;C=AA.c;D=AA.d;else;A=AA;end
end

DX1=D;D=0*D;
[A1,B1,C1,D1,A2,B2,C2,D2,M] = stabproj(A,B,C,D);
A2=A2+tol*eye(size(A2));

LC=lyap(-A2,B2*B2');LO=lyap(-A2',C2'*C2);
N1=LO*LC; sigr=sqrt(max(abs(eig(N1))));
if sigr<tol
   lam=0;
   AX=A1;BX=B1;CX=C1;DX=D1;
else
   lam=1/sigr*(1-tol/10);
   B2=B2*sqrt(lam);C2=C2*sqrt(lam);
   LO=LO*lam;LC=LC*lam;N1=N1*lam*lam;
   N=inv(eye(size(N1))-N1);
   a4=-A2';b4=N*LO*B2;c4=B2';d4=eye(size(c4*b4));
   a3=-A2';b3=N*C2';c3=-B2';d3=zeros(size(c3*b3));
   a2=A2;b2=N'*B2;c2=C2;d2=zeros(size(c2*b2));
   a1=A2;b1=-LC*N*C2';c1=C2;d1=eye(size(c1*b1));
   
   if Y==0
      [ai,bi,ci,di]=ssinv(a4,b4,c4,d4);
      ak=a2;bk=b2;ck=c2;dk=d2;
   else
      [ai,bi,ci,di]=addss(a4,b4,c4,d4,a3,b3*Y,c3,d3*Y);
      [ai,bi,ci,di]=ssinv(ai,bi,ci,di);
      [ak,bk,ck,dk]=addss(a2,b2,c2,d2,a1,b1*Y,c1,d1*Y);
   end
   
   
   [a,b,c,d]=series(ai,bi,ci,di,ak,bk,ck,dk);
   [AX,BX,CX,DX]=addss(A2,B2,C2,D2,a,-b,c,-d);
   [AX,BX,CX,DX]=addss(A1,B1,C1,D1,AX,BX/sqrt(lam),CX/sqrt(lam),DX);
end

[AX,BX,CX,DX]=stabproj(AX,BX,CX,DX);
DX=DX1;
eA=abs(eig(AX));cut=length(find(eA<0.01/tol));
if cut==0 & nosf == 0
  DX=DX-CX*inv(AX)*BX;
  AX=[];BX=[];CX=[];
else
  if cut<length(eA) & nosf == 0
    [AX,BX,CX,DX,AY,BY,CY,DY]=slowfast(AX,BX,CX,DX,cut);
    DCY=DY-CY*inv(AY)*BY;
    DX=DX+DCY;
  end
end

eA=abs(eig(AX));cut=length(find(eA<100*tol));
if cut==length(eA) & nosf == 0
    DX=DX;
    AX=[];BX=[];CX=[];
else
    if cut>0 & nosf == 0
        [AY,BY,CY,DY,AX,BX,CX,DX]=slowfast(AX,BX,CX,DX,cut);
        DX=DX+DY;
    end
end

temp=version;
if str2num(temp(1))>=6;
    if nargout <=2; AX=ss(AX,BX,CX,DX);BX=sigr;end
end
