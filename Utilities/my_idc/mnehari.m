function [AX,BX,CX,DX,sigr]=mnehari(G,B,C,D,tol);
% [AX,BX,CX,DX,sigr]=mnehari(A,B,C,D,tol);
% Solves the Nehari approximation problem
%      min_X ||R - X||_inf
%        s.t., R in L_inf, X in H_inf
%   R is in [A,B,C,D]
%   X is in [AX,BX,CX,DX]; sigr is the minimum distance
%   tol = Hankel inversion tolerance; high frequency modes are
%         cut beyond 0.1/tol. To bypass this reduction, give tol as
%         a negative number

nosf=0;
if nargin < 5; tol=1.e-6;        end
if tol==0;     tol=1e-6;         end
if tol<0;      tol=-tol; nosf=1; end

if isobject(G);
    G=ss(G);A=G.a;B=G.b;C=G.c;D=G.d;
else
    A=G;
end
    
[nout,nin]=size(D);
if nout > nin
    B=[B,0*C(1:nout-nin,:)'];D=[D,0*D(:,1:nout-nin)];
elseif nin > nout
    C=[C;0*B(:,1:nin-nout)'];D=[D;0*D(1:nin-nout,:)];
end
DX1=D;D=0*D;

G=ss(A,B,C,D);
[Gs,Gu]=stabproj(G); Gs=Gs+DX1;

Gu=Gu'; A2 = Gu.a; B2 = Gu.b; C2 = Gu.c; D2 = Gu.d;
%[A2,B2,C2]=obalreal(A2,B2,C2);
P=lyap(A2,B2*B2');Q=lyap(A2',C2'*C2);
QP=Q*P; sigr=sqrt(max(abs(eig(QP))));
if sigr<tol
    Y=ss([],[],[],0);
else
    m = (sigr*(1+tol))^2;
    GM = QP - m*eye(size(QP));
    AH = GM\(m*A2'+Q*A2*P);
    BH = GM\(Q*B2);
    CH = -C2*P;
    DH = -D2;
    YS = ss(AH,BH,CH,DH);
    Y = -YS';
end

eA=abs(eig(Y.a));cut=length(find(eA<0.1/tol));
if cut==0 & nosf == 0
    Y=Y.d-(Y.c)*inv(Y.a)*Y.b;
    Y=ss([],[],[],Y);
else
    if cut<length(eA) & nosf == 0
        [Ys,Yf]=slowfast(Y,cut);
        Df=Yf.d-Yf.c*inv(Yf.a)*Yf.b;
        Y=Ys+Df;
    end
end

eA=abs(eig(Y.a));cut=length(find(eA<10*tol));
if cut==length(eA) & nosf == 0
    Y=ss([],[],[],Y.d);
else
    if cut>0 & nosf ==0
        [Ys,Yf]=slowfast(Y,cut);
        Ds=Ys.d;
        Y=Yf+Ds;
    end
end

GX = minreal(Gs+Y);
GX = GX(1:nout,1:nin);
AX = GX.a; BX = GX.b; CX = GX.c; DX = GX.d;

if nargout <3
    BX=sigr;
    AX=GX;
end
