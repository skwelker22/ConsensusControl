function [f_sav,sc]=mvid(u,y,F_bw,F_ord,bw_cl0,P_num,P_den,IC_pts,m_red,fl_nr,fl_ns,LSMN0,SCALES0,FLAGS0,U_SCAL0,Y_SCAL0,targ_T,targ_S);

%function [file_ns,succ]=mvid(u,y,F_bw,F_ord,bw_cl,P_num,P_den,
%  IC_pts,m_red,fl_nr,fl_ns,LSMN,SCALES,FLAGS,U_SCAL,Y_SCAL,targ_T,targ_S);
% function MVID for the parametric identification of a MIMO
%  response contained in [u,y]
% Inputs:
%   F_ord = filter order: order of filters used for each MISO id
%   F_bw = filter BW: bandwidth of filters (i.e., 1/(s+BW)^ord).
%           BW*time-step should be less than 2 to 5. Typical values
%           for BW are around the expected closed-loop crossover 
%           frequency.
%   targ_T/S define the expected closed-loop bandwidth and
%      roll-off rates (in equiv. poles/zeros) of the closed loop
%      sensitivities 
%   file_ns = file_name to save ID data and parameters 
%            (2 or more chars)
%   in LSMN: error threshold and reduction weights: (min-dist. option)
%      e.g., [e,a,b] ([1,10,5]) searches within parameters 
%      that produce error <= (1+e)*LS_err, to minimize 
%      off-diagonal numerators and difference of denominator
%      coefficients from the filter, with weights 10 and 5.
%      Effectively attempts to diagonalize/stabilize the id-estimate.
%      typical selection: cutoff = [.5,3,3].
%   P_num,P_den:
%      prefilter poles and zeros: applied to both u and y to emphasize
%      the fit in a frequency region; they are entered as a vector
%      of corner frequencies; e.g., zeros=[.1 10], poles=[.4,4]
%      filters the data by (s+.1)(s+10)/(s+.4)(s+4); equal scalar
%      entries (e.g., 1 & 1) bypass the prefiltering step. 
%   IC_pts: number of points to be averaged for setting the prefilter
%      initial conditions.
%      Errors and uncertainties are computed from the unweighted u,y.
%   in SCALES: sampling (uniform), and time/input scaling factors.
%      The last two are transparent but may affect the location of 
%      the identified poles and zeros; their use is to find models
%      with similar error properties but easier to control pole-zero 
%      structure (e.g., try time scale -dilation- of 10). 
%    During manual Model Reduction (flag>=2): 
%      Exponential shifts may be used to stabilize the
%      model and allow the reduction of artificial RHP cancellations;
%      typically, use 0.1+max[real[pole]] where the latter is 
%      displayed immediately before the input request. 
%

%   KST-2/98-v.2

Y_SCALu=[];U_SCALu=[];Atot=[];
frun=[];AUNC=[];MUNC=[];W_ERR=[];t_T=[];t_S=[];sc=[];
if FLAGS0(1)==0;eval(['load ',fl_nr]);end
FLAGS=FLAGS0;SCALES=SCALES0;LSMN=LSMN0;bw_cl=bw_cl0;U_SCAL=U_SCAL0;Y_SCAL=Y_SCAL0;
hold off, format short e,
[nn,ninp]=size(u);[nn,noutp]=size(y);
valid=FLAGS(1);flag=FLAGS(2);
if length(FLAGS)>=3;ZER_MEAN=FLAGS(3);else;ZER_MEAN=0;end
stp=SCALES(1);tu_scal=SCALES(2:3);
fpts=SCALES(4);Mwin=SCALES(5);
nc=LSMN(1);red_meth=LSMN(2);
cutoff=LSMN(3:length(LSMN));
if length(cutoff)<5;cutoff(5)=1;end
if ZER_MEAN ==1
    nan_ind=[find(isnan(u(:,1)));nn+1];
    n_nan=length(nan_ind);f1_1=1;
    for i_nan = 1:n_nan
        ind_nan=[f1_1:nan_ind(i_nan)-1];
        f1_1=nan_ind(i_nan)+1;
        O_V = ones(length(u(ind_nan,:)),1);
        u(ind_nan,:)=u(ind_nan,:)-O_V*mean(u(ind_nan,:));
        y(ind_nan,:)=y(ind_nan,:)-O_V*mean(y(ind_nan,:));
    end
end

if valid==0 
    scalinp=1;t_dil=1;
else
    scalinp=tu_scal(2);t_dil=tu_scal(1);stp=stp/t_dil;
    F_bw=F_bw*t_dil;P_num=P_num*t_dil;P_den=P_den*t_dil;
end
while length(F_ord)<noutp,
    F_ord=[F_ord,F_ord(length(F_ord))];
end
while length(F_bw)<noutp,
    F_bw=[F_bw,F_bw(length(F_bw))];
end
fil_BW=log10(min(F_bw));
frw=logspace(-2+fil_BW,2+fil_BW,75)';
if length(cutoff)<5;
    while length(cutoff)<4;cutoff=[cutoff,1];end
end
if (cutoff(1))>0;cutoff(1)=max(cutoff(1),1.e-4);end
if (cutoff(2))<1.e-3;cutoff(2)=1.e-3;end
if (cutoff(3))<1.e-3;cutoff(3)=1.e-3;end
if (cutoff(4))<1.e-3;cutoff(4)=1.e-3;end
if (cutoff(5))<1.e-3;cutoff(5)=1.e-3;end

%------------------------- 
t=[0:nn-1]'*stp;clf;
subplot(121),plot(t,y(:,1));title('1st output')
subplot(122),plot(t,u(:,1));title('1st input')
pause(2)
%--------q2y=1;del0=.25;norflag=0;r1=1;r2=1;

%-------------------------
if valid ~=0
    
    %---------------------------------------- filters 
    u=u*diag(U_SCAL);y=y*diag(Y_SCAL);
    nuf=poly([-P_num])*P_den(1)/P_num(1); def=poly([-P_den]);
    uf=u;yf=y;
    if P_num ~= P_den
        uf=mimofilt(nuf,def,u,t,IC_pts);
        yf=mimofilt(nuf,def,y,t,IC_pts);
    end
end
subplot(121),plot(t,yf(:,1));title('1st scaled output')
subplot(122),plot(t,uf(:,1));title('1st scaled input')

%--------------------------------------------- balancing
if valid ~=0
    CCC=[];FFF=[];num=[10];
    for i=1:noutp
        den=poly(-ones(1,F_ord(i))*F_bw(i));n=F_ord(i);
        [f,q,cc,dd]=tf2ss(num,den);
        W=lyap(f,q*q');[U,S,V]=svd(W);TR=sqrt(inv(S))*U';
        f=TR*f*inv(TR);q=TR*q;
        CCC=[CCC,zeros(i-1,n);zeros(1,length(FFF)),q'];
        FFF=[FFF,zeros(length(FFF),n);zeros(n,length(FFF)),f'];
    end
    %----------------------------------------------------------
    disp('----------PARAMETER ESTIMATION----------------------')
    %----------------------------------------------------------
    if valid ==2
        [thx,err,nan_ind]=...
            par_estx(ninp,noutp,F_ord,FFF,CCC,uf,yf,t,cutoff,nc,red_meth);
        [Atot,Btot,Ctot,Dtot,Ftot,Bytot,Dytot,X0]=...
            m_xtrx(FFF,CCC,thx,ninp,noutp,F_ord,nan_ind);
    elseif valid==3
        %      par_esta
        [thx,err,nan_ind]=...
            par_esta(ninp,noutp,F_ord,FFF,CCC,uf,yf,t,cutoff,nc,red_meth);
        [Atot,Btot,Ctot,Dtot,Ftot,Bytot,Dytot,X0]=...
            mod_xtra(FFF,CCC,thx,ninp,noutp,F_ord,nan_ind);
    else
        [thx,err,nan_ind]=...
            par_est(ninp,noutp,F_ord,FFF,CCC,uf,yf,t,cutoff,nc,red_meth);
        [Atot,Btot,Ctot,Dtot,Ftot,Bytot,Dytot,X0]=...
            m_xtr(FFF,CCC,thx,ninp,noutp,F_ord,nan_ind);
    end
    clf;plot(t,err);title('Estimation Errors');pause(2)
    
    if flag>=0
        disp('Identified plant eigenvalues (scaled)')
        EAT=eig(Atot)
        disp('Identified plant transmission zeros (scaled)')
        TAT=tzero(Atot,Btot,Ctot,Dtot)
    end
    %----------------------------------------------------------
    if m_red(1)~=0
        disp('----------MODEL REDUCTION----------------------')
        %----------------------------------------------------------
        [Apr,Bprs,Cprs,Dprs,Fpr,Byprs,Dyprs,frw,MR]=...
            id_m_red(Atot,Btot,Ctot,Dtot,Ftot,Bytot,Dytot,F_bw,flag,m_red);
    else
        Apr=Atot;Bprs=Btot;Cprs=Ctot;Dprs=Dtot;Fpr=Ftot;Byprs=Bytot;Dyprs=Dytot;
        MR=[];
    end
    %----------- recover i/o scaling
    Bpr=Bprs*(diag(U_SCAL));Cpr=inv(diag(Y_SCAL))*Cprs;
    Dpr=inv(diag(Y_SCAL))*Dprs*(diag(U_SCAL));
    Bypr=Byprs*(diag(Y_SCAL));  Dypr=inv(diag(Y_SCAL))*Dyprs*(diag(Y_SCAL));
    u=u*inv(diag(U_SCAL));y=y*inv(diag(Y_SCAL));
end
clf
%----------------------------------

sspla=sv3_5(Apr,Bpr,Cpr,Dpr,1,frw);
loglog(frw,sspla);title('Identified Plant Singular Values');
if flag>1;pause;else;pause(2);end
%----------------------------------------------------------
disp('----------IC ESTIMATION (reduced)----------------------')
%----------------------------------------------------------

X0i=X0;
X0=[];     i_nan1=1;
for i_nan=1:length(nan_ind);
    ind_nan=[i_nan1:nan_ind(i_nan)-1];
    i_nan1=nan_ind(i_nan)+1;
    UY_nan=[u(ind_nan,:),y(ind_nan,:)];
    Y_nan=y(ind_nan,:);T_nan=t(ind_nan);T_nan=T_nan-T_nan(1);
    xx0=icestim(Fpr,[Bpr-Bypr*Dpr,Bypr],Cpr,[Dpr,Dypr],UY_nan,Y_nan,T_nan);
    X0=[X0,xx0];
end
if length(Atot)==length(Apr)
    disp(['IC normalized difference = ',num2str(norm(X0-X0i)/norm(X0))]);  
    disp(['IC norms = ',num2str(norm(X0i)),'  ',num2str(norm(X0))]);  
    if flag >0
        quer_ic = input('Select IC: 0=param-est., 1=post-est. [1]  ');
        if length(quer_ic)<1;quer_ic=1;end
    else
        quer_ic=1;
    end
    if quer_ic == 0;X0=X0i;end
end

if length(Apr)~=length(Atot) & flag >=0
    disp('Identified reduced plant eigenvalues (scaled)')
    EAT=eig(Apr)
    disp('Identified reduced plant transmission zeros (scaled)')
    TAT=tzero(Apr,Bpr,Cpr,Dpr)
end

disp('Identified DC gain')
DC=Dpr-Cpr*inv(Apr)*Bpr
%     disp('Identified scaled DC gain')
%     (diag(Y_SCAL))*DC*inv(diag(U_SCAL))
%----------------------------------------------------------
disp('----------IDENTIFICATION TESTS----------------------')
%----------------------------------------------------------
if flag>=1
    querytem=input('Compare output time responses? (0=no) [1]  ');
    if isempty(querytem);querytem=1;end
else
    querytem=1;
end
if querytem ~=0
    Y=lsim_nan(Apr,Bpr,Cpr,Dpr,u,t,X0);
    for i=1:noutp
        clf,hold off,subplot(111)
        plot(t,y(:,i),'r',t,Y(:,i),'g');
        title('Actual (r) vs. Identified (g) system comparison');
        if flag>=2;pause;else;pause(2);end
    end
    plot(t,y-Y);title('Additive Identification Error');
    if flag>=2;pause;else;pause(2);end
end
%----------------------------------------------------------
disp('----------UNCERTAINTY ESTIMATION----------------------')
%----------------------------------------------------------
frun=[];AUNC=[];MUNC=[];succ=1;
Bprs=Bpr*inv(diag(U_SCAL));Cprs=diag(Y_SCAL)*Cpr;
Dprs=diag(Y_SCAL)*Dpr*inv(diag(U_SCAL));
Byprs=Bypr*inv(diag(Y_SCAL));Dyprs=(diag(Y_SCAL))*Dypr*inv(diag(Y_SCAL));
if flag>=1
    querytem=input('Perform uncertainty analysis? (0=no) [1]  ');
    if isempty(querytem);querytem=1;end
else
    querytem=1;
end
if querytem ~=0
    i_nnan=find(~isnan(u(:,1)));
    quer=input('Compute new i/o scales? (1=err,-1=out) [1,-1]  ');
    if length(quer)<1;quer=1;end
    if length(quer)<2;
        if quer~=0;quer=[quer,-1];else;quer=[0,0];end
    end
    u_scal=ones(ninp,1);y_scal=ones(noutp,1);
    while norm(quer)~=0
        u_scal=ones(ninp,1);y_scal=ones(noutp,1);
        if quer(1)~=0
            u_rms=diag((u(i_nnan,:)*diag(U_SCAL))'*(u(i_nnan,:)*diag(U_SCAL)));
            min_u=min(u_rms); u_scal=sqrt(min_u./u_rms);
        end
        if quer(2) == 1
            wu=lsim_nan(Fpr,Bprs-Byprs*Dprs,Cprs,Dprs,u*diag(U_SCAL),t,X0);
            wy=lsim_nan(Fpr,Byprs,Cprs,Dyprs,y*diag(Y_SCAL),t);
            Ferr=y*diag(Y_SCAL)-(wu+wy);
            err_rms=diag(Ferr(i_nnan,:)'*Ferr(i_nnan,:));
            av_err=mean(err_rms); y_scal=sqrt(av_err./err_rms);
        elseif quer(2) == -1
            y_rms=diag((y(i_nnan,:)*diag(Y_SCAL))'*(y(i_nnan,:)*diag(Y_SCAL)));
            min_y=min(y_rms);y_scal=sqrt(min_y./y_rms);
        end
        disp('Input/Output scale adjustments')
        disp(u_scal');disp(y_scal')
        quer=input('Compute new i/o scales? (1=err,-1=out) [0,0]  ');
        if length(quer)<1;quer=0;end
        if length(quer)<2;quer=[quer,0];end
        if quer(1)==-1;quer=[0,0];u_scal=ones(ninp,1);y_scal=ones(noutp,1);end
        
    end
    U_SCALu=U_SCAL.*u_scal;Y_SCALu=Y_SCAL.*y_scal;
    Bprs=Bpr*inv(diag(U_SCALu));Cprs=diag(Y_SCALu)*Cpr;
    Dprs=diag(Y_SCALu)*Dpr*inv(diag(U_SCALu));
    Byprs=Bypr*inv(diag(Y_SCALu));
    Dyprs=(diag(Y_SCALu))*Dypr*inv(diag(Y_SCALu));
    us=u*diag(U_SCALu);ys=y*diag(Y_SCALu);
    
    %    Dyi=inv(eye(size(Dyprs))-Dyprs);
    %%    Bprs=Bprs+Byprs*Dyi*Dprs;
    %    Cprs=Dyi*Cprs;
    %    Dprs=Dyi*Dprs;
    
    if flag>=1
        plot(t,us);title('scaled input');pause(1)
        plot(t,ys);title('scaled output');pause(1)
        sspsc=sv3_5(Apr,Bprs,Cprs,Dprs,1,frw);
        loglog(frw,sspsc);title('Scaled Plant Singular Values');
        pause(1)
    end
    
    %------------------------------- recover time dilation
    Apr=Apr/t_dil;Bpr=Bpr/t_dil;Fpr=Fpr/t_dil;Bypr=Bypr/t_dil;
    Bprs=Bprs/t_dil;Byprs=Byprs/t_dil;
    stp=stp*t_dil;frw=frw/t_dil;t=t*t_dil;DC=Dpr-Cpr*inv(Apr)*Bpr;
    F_bw=F_bw/t_dil;P_num=P_num/t_dil;P_den=P_den/t_dil;
    
    % ----------------     
    [frun,AUNC,MUNC,W_ERR,t_T,t_S,sc]=uncertn(Apr,Fpr,Bprs,Cprs,Dprs,Byprs,Dyprs,X0,us,ys,t,fpts,flag);
end
%------------------------------- save results

FLAGS=[valid,flag,ZER_MEAN];tu_scal=[t_dil,scalinp];
SCALES=[stp,tu_scal,fpts,Mwin];
LSMN=[nc,red_meth,cutoff];

if flag >=2 | length(fl_ns)<2
    savsys=input(['Give filename to save ID results; current: ',fl_ns,'  '],'s');
else
    savsys=fl_ns;
end

if length(savsys)>=2
    eval(['save ',savsys,' Apr Bpr Cpr Dpr Fpr Bypr Dypr X0 DC frw sspla MR frun MUNC AUNC F_bw F_ord P_num P_den IC_pts m_red LSMN SCALES FLAGS U_SCAL Y_SCAL U_SCALu Y_SCALu t_T t_S W_ERR bw_cl'])
    f_sav=savsys;
else
    f_sav=[];
end
