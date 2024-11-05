function [AX,BX,CX,DX,sigr]=nehari(A,B,C,D,tol);
% [AX,BX,CX,DX,sigr]=nehari(A,B,C,D,tol);
% Solves the Nehari approximation problem
%      min_X ||R - X||_inf
%        s.t., R in L_inf, X in H_inf 
%   R is in [A,B,C,D]
%   X is in [AX,BX,CX,DX]; sigr is the minimum distance


if nargin < 5;tol=1.e-6;end

[A1,B1,C1,D1,A2,B2,C2,D2,M] = STABPROJ(A,B,C,D);
D1=D1+D2;D2=0*D2;A2=A2+tol*eye(size(A2));
[Ab,Bb,Cb,G,T] = OBALREAL(-A2',C2',B2');
sigr=max(G);A2=-Ab';B2=Cb';C2=Bb';
sig=1/sigr*(1-(tol)/100);
NS=inv(eye(length(G),length(G))-sig*sig*diag(G)*diag(G));
AE=[A2,sig*NS'*B2*B2';0*A2,-(A2'+sig*sig*NS*diag(G)*B2*B2')];
BE=[sqrt(sig)*NS'*B2;-sqrt(sig^3)*NS*diag(G)*B2];
CE=[C2,0*C2];DE=0*CE*BE;
[AX,BX,CX,DX]=addss(A2,sig*B2,C2,sig*D2,AE,BE,-CE,-DE);
[AX,BX,CX,DX]=addss(A1,B1,C1,D1,AX,BX,CX,DX);
[AX,BX,CX,DX,A2,B2,C2,D2,M] = STABPROJ(AX,BX,CX,DX);
DX=DX+D2;

