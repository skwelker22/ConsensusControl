function [S_NL,S_DL,S_NR,S_DR,S_XL,S_YL,S_XR,S_YR]=doub_cop(A,B,C,D,shif,gain);

%function [S_NL,S_DL,S_NR,S_DR,S_XL,S_YL,S_XR,S_YR]=...
%                doub_cop(A,B,C,D,shif,gain);
% Find doubly coprime factorizations for a system [A,B,C,D]
% P = NR DRinv = DLinv NL
% [ XL  -YL ] [ DR   YR ]
% [-NL   DL ] [ NR   XR ] = I
%    Unpack using branch: [as,bs,cs,ds]=branch(SS_s);
%    shif (optional) contains the guaranteed stability margin of
%        the coprime factors
%    gain (optional) contains the LQ control weight used in the stabilization

% S_ = mksys(A,B,C,D,'ss');

if nargin <5; shif=1;end
if nargin <6; gain=1;end

nA=length(A);[nout,ninp]=size(D);
IA=eye(nA,nA);IO=gain*eye(nout,nout);OI=gain*eye(ninp,ninp);
     Lpr=lqr(A'+shif*IA,C',B*B',IO);Lpr=Lpr';
     Kpr=lqr(A+shif*IA,B,C'*C,OI);

a=A-Lpr*C;c=C;b=B-Lpr*D;d=D;
S_NL = mksys(a,b,c,d,'ss');

a=A-Lpr*C;c=C;b=-Lpr;d=eye(nout,nout);
S_DL = mksys(a,b,c,d,'ss');

a=A-B*Kpr;c=C-D*Kpr;b=B;d=D;
S_NR = mksys(a,b,c,d,'ss');

a=A-B*Kpr;c=-Kpr;b=B;d=eye(ninp,ninp);
S_DR = mksys(a,b,c,d,'ss');

a=A-Lpr*C;c=-Kpr;b=-(B-Lpr*D);d=eye(size(c*b));
S_XL = mksys(a,b,c,d,'ss');

a=A-Lpr*C;c=-Kpr;b=Lpr;d=0*(c*b);
S_YL = mksys(a,b,c,d,'ss');

a=A-B*Kpr;c=-Kpr;b=Lpr;d=0*c*b;
S_YR = mksys(a,b,c,d,'ss');

a=A-B*Kpr;c=C-D*Kpr;b=Lpr;d=eye(size(c*b));
S_XR = mksys(a,b,c,d,'ss');

