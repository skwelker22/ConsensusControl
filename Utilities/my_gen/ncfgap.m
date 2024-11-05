function [NS,MS,Q,sigr]=ncfgap(N1,M1,N2,M2,SW_N,red_tol,qslow,mrtol);

%function [D_S,N_S,Q,sigr]=ncfgap(N_1,M_1,N_2,M_2,SW_N,red_tol,qslow);
% find the directional gap distance between Sys1 and Sys2;
% also returns the optimal Q minimizing G1 - Q G2

if nargin < 5, SW_N = []; end
if nargin < 6, red_tol = 1e-6; end
if nargin < 7, qslow = 1e-4; end
if nargin < 8, mrtol = 1e-8; end

if isempty(M1), [N1,M1]=ncf(N1); end
if isempty(M2), [N2,M2]=ncf(N2); end

if ~isempty(SW_N), end
%G1 = [N1 , M1]; G2 = [N2 , M2];
if ~isempty(SW_N), end

G1 = [-M1 , N1]; G2 = [-M2 , N2];
P2=minreal(inv(M2)*N2,mrtol); G2r=ncf(P2,1);
P1=minreal(inv(M1)*N1,mrtol); G1r=ncf(P1,1);
Gd = G1*G2r; Gr = G1*G2';

gmax=norm(G1,inf); gmin=norm(Gd,inf);    [gmax,gmin]
while gmax-gmin>0.0001
    gam=(gmax+gmin)/2; G0=(sfr(minreal(Gd,mrtol)/gam))*gam;
    R=minreal(Gr*inv(G0),mrtol);[Qn,gtest]=nehari(R);
    if gtest < 1; gmax=gam; else gmin=gam;end
end
Q=minreal(Qn*G0,mrtol);
cut=length(find(abs(eig(Q))<qslow));[qs,qf]=slowfast(Q-Q.d,cut);qf=qf+Q.d;
Q=h_sysred(qf,[],[],0,0,[3,red_tol]);
NS = minreal(N1-Q*N2,mrtol); MS = minreal(M1-Q*M2,mrtol);
NS=h_sysred(NS,[],[],0,0,[3,red_tol]);MS=h_sysred(MS,[],[],0,0,[3,red_tol]);

nugap = norm(G2*G1r,inf);

%n=hankelsv(G1); n=max(n);
Wc=gram(G1,'c');Wo=gram(G1,'o');
n=max(eig(Wc*Wo));
optgap1=sqrt(1-n);
%max_unc(P1)
Wc=gram(G2,'c');Wo=gram(G2,'o');
n=max(eig(Wc*Wo));
optgap2=sqrt(1-n);

sigr = [gmax nugap optgap1 optgap2];

if nargout == 0;NS=sigr;disp(' ');disp('gmax,   nugap,   max_unc1,   maxunc2');end
