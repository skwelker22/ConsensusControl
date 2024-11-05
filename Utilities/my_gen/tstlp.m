function [frun,ssu,ssy]=tstlp(t,u,y);
% Usage function [ssu,ssy]=fcesens(pn_,pd_,c_,frw);
%   computes the singular values of the closed loop sensitivities
%   from the coprime factor error e to the plant input and output 
%   (u,y).  pn_,pd_,c_ are the mksys-packed coprime factors and 
%   controller
%                 u          |e             y
%          --> c_ --> pn_ -->o--> inv(pd_) -->
%          |                                 |
%          |_________________________________|
%

inid=input('Inner loop id  ','s');
incon=input('Inner loop controller  ','s');

eval(['load ' inid]);
eval(['load ' incon]);

[AN,BN,CN,DN,AD,BD,CD,DD]=norm_cop(Apr,Bpr,Cpr,Dpr);
c_=mksys(Acr,Bcr,Ccr,Dcr,'ss');
pn_=mksys(AN,BN,CN,DN,'ss');
pd_=mksys(AD,BD,CD,DD,'ss');

[ssu,ssy]=cfesens(pn_,pd_,c_,frun);

[noutp,ninp]=size(Dpr);
[FN,CN,BUN,DUN,BYN,DYN]=normod(Apr,Bpr,Cpr,Dpr);
wu=lsim(FN,[BUN,BYN],CN,[DUN,DYN],[u,y],t);
Ferr=[];AAA=[];
  for i=1:noutp
    Ferr=[Ferr;y(:,i)-(wu(:,i))];
    w0=lsim(FN',zeros(length(FN),1),eye(length(FN)),...
            zeros(length(FN),1),0*t,t,CN(i,:)');
    AAA=[AAA;w0];
  end
X0=inv(AAA'*AAA)*AAA'*Ferr;


% construct errors
wu=lsim(FN,[BUN,BYN],CN,[DUN,DYN],[u,y],t,X0);
Ferr=y-(wu);
clg
plot(t,Ferr),title('Estimation Error'),pause

disp('Uncertainty bound estimation ');
fpoints=80;stp=t(2)-t(1);

[frunz,UNCN]=unc_bndf(Ferr,u,FN,BYN,CN,DYN,0,fpoints,stp);
if norm(frun-frunz)>1.e-10;disp('Error in the frequency range');end
[frunz,UNCD]=unc_bndf(Ferr,y,FN,BYN,CN,DYN,0,fpoints,stp);
[frunz,UNCND]=unc_bndf(Ferr,u,y,BYN,CN,DYN,-1,fpoints,stp);

UNCNi=1.0./UNCN;UNCDi=1.0./UNCD;
ssum=max(ssu)';ssym=max(ssy)';
loglog(frun,UNCNi,'r',frun,ssu);
grid;title('Input Based Uncertainty')
pause
loglog(frun,UNCDi,'r',frun,ssy);
grid;title('Output based Uncertainty')
pause
loglog(frun,ssum.*UNCN,'y',frun,ssym.*UNCD,'r')
grid;title('SGT Test')
pause
loglog(frun,sqrt((ssum.^2+ssym.^2).*(UNCND(:,1).^2+UNCND(:,2).^2)))
grid;title('SGT Test II')
pause
x='done'

