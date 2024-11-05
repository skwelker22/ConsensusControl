% Script file MIMOID2 for the parametric identification of a MIMO
%  response contained in [u,y]
% Inputs:
%    filter order: order of filters used for each MISO id
%    filter BW: bandwidth of filters (i.e., 1/(s+BW)^ord)
%      This version iterates the aux.filters to minimize
%      the low-freq. additive uncertainty
%    error threshold and reduction weights: (min-dist. option)
%      e.g., [e,a,b] ([1,10,5]) searches within parameters 
%      that produce error <= (1+e)*LS_err, to minimize 
%      off-diagonal numerators and difference of denominator
%      coefficients from the filter, with weights 10 and 5.
%      Effectively attempts to diagonalize/stabilize the id-estimate.
%    prefilter poles and zeros: applied to both u and y; errors
%      and uncertainties are computed from the unweighted u,y.
%

format short e,Mwindow=4;scalinp=1;
stp0=input('time step [0.5/60]  ');
   [nn,ninp]=size(u);[nn1,noutp]=size(y);
     if nn ~= nn1;
     disp('id error: lengths of y and u must be equal ');
     return;end
   if length(stp0) > 0, 
     if stp0~=0, stp=stp0;end
   else, stp=.5/60
   end
hold off, 
t=[0:nn-1]'*stp;clg;subplot(121),plot(t,y(:,1));title('1st output')
subplot(122),plot(t,u(:,1));title('1st input')
q2y=1;del0=.25;norflag=0;r1=1;r2=1;fpoints=80;
valid=input('identification (>0) or validation (0) [1]  ');
if length(valid)<1,valid=1;end
def_flag=input('use default selections? (1=yes) [0]  ');
if length(def_flag)<1,def_flag=0;end

if valid ~=0
  filord=input('filter order  [3] ');
    if length(filord)<1,filord=3;end
    while length(filord)<noutp,
       filord=[filord,filord(length(filord))];
    end
  filpol=input('filter bandwidth [5]  ');
    if length(filpol)<1,filpol=5;end
    while length(filpol)<noutp,
       filpol=[filpol,filpol(length(filpol))];
    end
  fil_BW=min(filpol);
    if def_flag == 0
    nc=input('number of constraints (0-biproper,1-str.proper) [0]  ');
      if length(nc)<1,nc=0;end
    cutoff=input('LS threshold and weights (stability, off-diag N/D,) [.5 1 1 1]  ');
      if length(cutoff)<1;cutoff=0.5;end
    cutoff=abs(cutoff);
      if length(cutoff)==1;cutoff=[cutoff,1,1,1];
      elseif length(cutoff)==2; cutoff=[cutoff,1,1];
      elseif length(cutoff)==3; cutoff=[cutoff,1];
      end
      if (cutoff(1))>0;cutoff(1)=max(cutoff(1),1.e-4);end
      if (cutoff(2))<1.e-3;cutoff(2)=1.e-3;end
      if (cutoff(3))<1.e-3;cutoff(3)=1.e-3;end
      if (cutoff(4))<1.e-3;cutoff(4)=1.e-3;end
    red_meth=input('DOF reduction method (1=SVD,2=min.dist.) [2]  ');
      if length(red_meth)<1;red_meth=2;end
    prenum=input('prefilter -zeros [1]  ');
    preden=input('prefilter -poles [1]  ');
    scalinp=input('input scaling factor [1]  ');
      if length(prenum)<1;prenum=1;end
      if length(preden)<1;preden=1;end
      if length(scalinp)<1;scalinp=1;end
  else
    nc=0;cutoff=[.5 1 1 1];prenum=1;preden=1;red_meth=2;scalinp=1;
  end


% filters 
    u=u*scalinp;
    nuf=poly([-prenum])*preden(1)/prenum(1); def=poly([-preden]);
    if prenum ~= preden
      uf=mimofilt(nuf,def,u,t);yf=mimofilt(nuf,def,y,t);
    else
      uf=u;yf=y;
    end
end


if valid ~=0
  CCC=[];FFF=[];num=[10];
    for i=1:noutp
%     balancing
      den=poly(-ones(1,filord(i))*filpol(i));n=filord(i);
      [f,q,cc,dd]=tf2ss(num,den);
      W=lyap(f,q*q');[U,S,V]=svd(W);TR=sqrt(inv(S))*U';
      f=TR*f*inv(TR);q=TR*q;
      CCC=[CCC,zeros(i-1,n);zeros(1,length(FFF)),q'];
      FFF=[FFF,zeros(length(FFF),n);zeros(n,length(FFF)),f'];
    end

  thx00=[];rdone=1;%y0=yf;

   while rdone ~= 0
     [thx,err,WWX]=par_est(ninp,noutp,filord,FFF,CCC,uf,yf,t,cutoff,nc,red_meth);
     [Atot,Btot,Ctot,Dtot,Ftot,Bytot,Dytot,X0]=...
                       mod_xtr(FFF,CCC,thx,ninp,noutp,filord);
     EAT=eig(Atot),TAT=tzero(Atot,Btot,Ctot,Dtot)
     querytem=input('Gain Feedback Evals ? (y~=0) ');
       if length(querytem) == 1 & ninp==noutp
         FAT=eig(Atot-Btot*Ctot*querytem)
       end

     [frun,AUNC]=unc_bndf(err,u,Ftot,Bytot,Ctot,Dytot,3,fpoints,stp);
     loglog(frun,AUNC);grid;title('Additive Uncertainty');pause
     SVN=sv3_5(Atot,Btot,Ctot,Dtot,1,frun);
       if length(thx00)>1, norm(thx-thx00),
          loglog(frun,SVN,'r',frun,SVO,'b');
          title('singular value comparison; red=new, blue=old');pause
          loglog(frun,AUNC,'r',frun,AUNCO,'b');grid;
          title('Additive Uncertainty comparison; red=new, blue=old');pause
       end
     thx00=thx;SVO=SVN;AUNCO=AUNC;

     querytem=input('Compare output time responses? (0=no) ');
       if querytem ~=0
         Y=lsim(Atot,Btot,Ctot,Dtot,u,t,X0);
         add_err=norm(y-Y)
           for i=1:noutp
             plot(t,y(:,i),'r',t,Y(:,i),'g');
             title('Actual (r) vs. Identified (g) system comparison');
             pause
           end
         plot(t,y-Y);title('Additive Identification Error');
         pause
       end

     rdone=input('Repeat ? (y~=0) ');
       if rdone ~=0
        [q1,q2,q3,q4,q5,q6,q7,q8,atem,btem,ctem,dtem]=...
           iofc(Ftot,Bytot,-Ctot,eye(size(Dytot))-Dytot);
           [aax,bbx,ccx,ddx]=ssinv(atem,btem,ctem,dtem);
           ss_dcor=sv3_5(aax,bbx,ccx,ddx,1,frun);
           loglog(frun,ss_dcor);grid;title('Adjustment svs');pause
 %       mfilpol=.01;
 %       Klq=lqr(Atot'+mfilpol/100*eye(size(Atot)),eye(size(Atot)),...
 %            eye(size(Atot)),10000*eye(size(Atot)));
 %        Ftot=Atot-Klq;
          FFF=aax;
       end
  end

  [Apr,Bpr,Cpr,Dpr,Fpr,Bypr,Dypr,frw,MR]=...
       mod_red(Atot,Btot,Ctot,Dtot,Ftot,Bytot,Dytot,filpol);

end
clg

u=u/scalinp;uf=uf/scalinp;Bpr=Bpr*scalinp;Dpr=Dpr*scalinp;MR=MR*scalinp;

sspla=sv3_5(Apr,Bpr,Cpr,Dpr,1,frw);
loglog(frw,sspla);title('Identified Plant Singular Values');pause

X0=icestim(Fpr,[Bpr-Bypr*Dpr,Bypr],Cpr,[Dpr,Dypr],[u,y],y,t);

disp('Identified reduced plant eigenvalues')
EAT=eig(Apr)
disp('Identified reduced plant transmission zeros')
TAT=tzero(Apr,Bpr,Cpr,Dpr)
  querytem=input('Gain Feedback Evals ? (y~=0) [0]  ');
     if length(querytem) == 1 & ninp==noutp
       FAT=eig(Apr-Bpr*Cpr*querytem)
     end

querytem=input('Compare output time responses? (0=no) ');
    if querytem ~=0
      Y=lsim(Apr,Bpr,Cpr,Dpr,u,t,X0);
        for i=1:noutp
          plot(t,y(:,i),'r',t,Y(:,i),'g');
          title('Actual (r) vs. Identified (g) system comparison');
          pause
        end
      plot(t,y-Y);title('Additive Identification Error');
      pause
    end

querytem=input('Perform frequency domain analysis ? (1=y) ');
     if querytem ==1
        FREQAN
     end

querytem=input('Perform uncertainty analysis ? (0=no) ');
   if querytem ~=0
[frun,AUNC,MUNC]=uncerf(Apr,Fpr,Bpr,Cpr,Dpr,Bypr,Dypr,X0,u,y,t,fpoints,stp,def_flag);
   end

disp('Other uncertainty computation options: UNCERTT,UNCERTL')
disp('If you want to save the id results, use:')
disp('save FID1xxx stp den cutoff Apr Bpr Cpr Dpr Fpr Bypr Dypr frw MR frun MUNC AUNC')
disp('save FID1fd frw MR fftfr Gest frg Gspa SDGspa') 

savsys=input('Give filename to save id results   ','s');
if length(savsys)>2
eval(['save ',savsys,' stp den cutoff Apr Bpr Cpr Dpr Fpr Bypr Dypr frw MR frun MUNC AUNC'])
end

