function [ssu,ssy]=cfesens(pn_,pd_,c_,frw);
% Usage function [ssu,ssy]=cfesens(pn_,pd_,c_,frw);
%   computes the singular values of the closed loop sensitivities
%   from the coprime factor error e to the plant input and output 
%   (u,y).  pn_,pd_,c_ are the mksys-packed coprime factors and 
%   controller
%                 u          |e             y
%          --> c_ --> pn_ -->o--> inv(pd_) -->
%          |                                 |
%          |_________________________________|
%

[An,Bn,Cn,Dn]=branch(pn_);
[Ad,Bd,Cd,Dd]=branch(pd_);
[Ac,Bc,Cc,Dc]=branch(c_);
[Ad,Bd,Cd,Dd]=ssinv(Ad,Bd,Cd,Dd);

nn=length(An);nd=length(Ad);nc=length(Ac);
[in,on]=size(Dn');[id,od]=size(Dd');[ic,oc]=size(Dc');

At=[Ac,zeros(nc,nd+nn);zeros(nn,nc),An,zeros(nn,nd);zeros(nd,nc+nn),Ad];
Bt=[Bc,zeros(nc,in+id);zeros(nn,ic),Bn,zeros(nn,id);zeros(nd,ic+in),Bd];
Ct=[Cc,zeros(oc,nd+nn);zeros(on,nc),Cn,zeros(on,nd);zeros(od,nc+nn),Cd];
Dt=[Dc,zeros(oc,id+in);zeros(on,ic),Dn,zeros(on,id);zeros(od,ic+in),Dd];

q1=[1:ic]';q2=[ic+1:ic+in]';q3=[ic+in+1:ic+in+id]';
q4=[1:oc]';q5=[oc+1:oc+on]';q6=[oc+on+1:oc+on+od]';
Q=[q1 q6;q2 q4;q3 q5];IN=q3;OUTu=[q4];OUTy=[q6];

[Au,Bu,Cu,Du] = CONNECT(At,Bt,Ct,Dt,Q,IN,OUTu);
[Ay,By,Cy,Dy] = CONNECT(At,Bt,Ct,Dt,Q,IN,OUTy);

ssu=sv3_5(Au,Bu,Cu,Du,1,frw);
ssy=sv3_5(Ay,By,Cy,Dy,1,frw);




