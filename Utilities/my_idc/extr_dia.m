function [Apn,Bpn,Cpn,Dpn,Noutp,Ninp]=extr_dia(Apr,Bpr,Cpr,Dpr,sch_tol);
%function [Apn,Bpn,Cpn,Dpn,Noutp,Ninp]=extr_dia(Apr,Bpr,Cpr,Dpr,sch_tol);
% Extraction of the diagonal of a system; 
% model reductions use w_sysred with tolerance sch_tol [def=1.e-6]

if nargin<5; sch_tol=1.e-6;end

[Noutp,Ninp]=size(Dpr);
if Noutp~=Ninp
  disp('** extr_dia warning: System not square: resetting I/O')
  NNX=min(Noutp,Ninp);Noutp=NNX;Ninp=NNX;
  Bpr=Bpr(:,1:NNX);Cpr=Cpr(1:NNX,:);Dpr=Dpr(1:NNX,1:NNX);
end

if Noutp >1
   IX=1;
   A_sys=mksys(Apr,Bpr(:,IX),Cpr(IX,:),Dpr(IX,IX),'ss');
   [Aprx,Bprx,Cprx,Dprx]=w_sysred(A_sys,[],[],0,[],[3,sch_tol]);
   IX=2;
   A_sys=mksys(Apr,Bpr(:,IX),Cpr(IX,:),Dpr(IX,IX),'ss');
   [Aprx2,Bprx2,Cprx2,Dprx2]=w_sysred(A_sys,[],[],0,[],[3,sch_tol]);
   [Apn,Bpn,Cpn,Dpn]=append(Aprx,Bprx,Cprx,Dprx,Aprx2,Bprx2,Cprx2,Dprx2);
else
   Apn=Apr;Bpn=Bpr;Cpn=Cpr;Dpn=Dpr;
end
if Noutp > 2
   for IX=3:Noutp
   A_sys=mksys(Apr,Bpr(:,IX),Cpr(IX,:),Dpr(IX,IX),'ss');
   [Aprx,Bprx,Cprx,Dprx]=w_sysred(A_sys,[],[],0,[],[3,sch_tol]);
   [Apn,Bpn,Cpn,Dpn]=append(Apn,Bpn,Cpn,Dpn,Aprx,Bprx,Cprx,Dprx);
   end
end

