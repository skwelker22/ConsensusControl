function [xx,S_cl,T_cl]=h_perfev(Acr,Bcr,Ccr,Dcr,Apr,Bpr,Cpr,Dpr,frw,flag,figno);
%function [xx,S_cl,T_cl]=h_perfev(Acr,Bcr,Ccr,Dcr,Apr,Bpr,Cpr,Dpr,frw,flag);
%----------- CLOSED LOOP EVALUATION -----------------

xx=1;if nargin<10;flag=1;end
if nargin<11;figno=1;end

[Alpgn,Blpgn,Clpgn,Dlpgn]=series(Acr,Bcr,Ccr,Dcr,Apr,Bpr,Cpr,Dpr);
[Acly,Bcly,Ccly,Dcly]=feedbk(Alpgn,Blpgn,Clpgn,Dlpgn,2);
[Aclu,Bclu,Cclu,Dclu]=feedbk(Acr,Bcr,Ccr,Dcr,3,Apr,Bpr,Cpr,Dpr);
[Acle,Bcle,Ccle,Dcle]=feedbk(Alpgn,Blpgn,Clpgn,Dlpgn,1);
[Apus,Bpus,Cpus,Dpus]=feedbk(Apr,Bpr,Cpr,Dpr,3,Acr,Bcr,Ccr,Dcr);
T_cl=mksys(Acly,Bcly,Ccly,Dcly,'ss');S_cl=mksys(Acle,Bcle,Ccle,Dcle,'ss');
if flag > 1;
    EAC=eig(Acr);
    disp(['max controller eigenvalue: ',num2str(max(real(EAC)))])
    TZAC=tzero(Acr,Bcr,Ccr,Dcr);
    disp(['max controller T-zeros: ',num2str(max(real(TZAC))),', ',num2str(min(abs(TZAC)))])
end
cllpeig=eig(Acle);
disp(['Max closed-loop poles:  ',num2str(max(real(cllpeig)))])
if max(real(cllpeig))>0
    disp('CONTROLLER REDUCTION ERROR:');
    disp('   REPEAT WITH A DIFFERENT METHOD AND/OR THRESHOLDS');
    xx=0;
end

itemp=1;
if itemp ~=0
    magcon=sv3_5(Acr,Bcr,Ccr,Dcr,1,frw);
figure(figno);clf
    loglog(frw,sat(magcon',1e5,1.e-5))
    title('Controller SV')
    magcl=sv3_5(Acly,Bcly,Ccly,Dcly,1,frw);
    magcle=sv3_5(Acle,Bcle,Ccle,Dcle,1,frw);
    %   magur=sv3_5(Aclu,Bclu,Cclu,Dclu,1,frw);
    max_T=max(max(magcl));max_S=max(max(magcle));
    disp(['max |T| = ',num2str(max_T),',   max |S| = ',num2str(max_S)])
figure(figno+1);clf
    loglog(frw,sat([magcl',magcle'],100,1.e-5))
    title('S&T sensitivities')
    bwT=frw(min(find(magcl(1,:)<.7)))
    if isempty(bwT);bwT=input('h_perfev: give closed loop bandwidth  ');
        if length(bwT)<1;bwT=1;end;end
    Timmax=50/bwT;
    if flag >=3;timx=input('give max sim time  ');else;timx=[];end
    if ~isempty(timx);Timmax=timx;end
    T=[0:.004*Timmax:Timmax]';nnt=length(T);[nnx,nny]=size(Dcly);
    uu=ones(nnt,nny);YY=zeros(nnt,nnx*nny);
    for i=1:nny
        YY(:,(i-1)*nnx+1:i*nnx)=step(Acly,Bcly,Ccly,Dcly,i,T);
    end
    Y=lsim(Acly,Bcly,Ccly,Dcly,uu,T);
figure(figno+2);clf
    plot(T,YY);title('step responses for each channel');
figure(figno+3);clf
    plot(T,Y);title('step response (all inputs)');
end
