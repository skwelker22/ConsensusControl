function [f_sav,succ,SYS]=mimoid(u,y,F_bw,F_ord,cl_spec0,P_num,P_den,IC_pts,m_red,fl_nr,fl_ns,LSMN0,SCALES0,FLAGS0,ssTi,ssDmi);

%function [f_sav,succ]=mimoid(u,y,F_bw,F_ord,cl_spec0,P_num,P_den,IC_pts,
%              m_red,fl_nr,fl_ns,LSMN0,SCALES0,FLAGS0,ssTi,ssDmi);
%
% function MIMOID for the parametric identification of a MIMO
%  response contained in [u,y]
% Inputs:
%   F_ord = filter order: order of filters used for each MISO id
%   F_bw = filter BW: bandwidth of filters (i.e., 1/(s+BW)^ord).
%           BW*time-step should be less than 2 to 5. Typical values
%           for BW are around the expected closed-loop crossover 
%           frequency.
%   cl_spec defines the expected closed-loop bandwidth and
%      roll-off rates (in equiv. poles/zeros) of the closed loop
%      sensitivities and Risk and Aggressiveness factors used
%      for auto-spec [T_bw,T_roll,S_bw,S_roll,RISK,AGRS].
%   file_ns = file_name to save ID data and parameters 
%            (2 or more chars)
%   in LSMN: error threshold and reduction weights: (min-dist. option)
%      e.g., [e,a,b] ([1,10,5]) searches within parameters 
%      that produce error <= (1+e)*LS_err, to minimize 
%      off-diagonal numerators and difference of denominator
%      coefficients from the filter, with weights 10 and 5.
%      Effectively attempts to diagonalize/stabilize the id-estimate.
%      typical selection: cutoff = [.1,3,3].
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
%   ssTi,ssDmi = inner loop complementary sv's and uncertainty (rows)
%
%   Called by the interactive function MIMOIDI which specifies the
%   data structures

%   KST-7/25/97-v.2

VarU=1;VarY=1; MR=[];

if nargin<15;ssDmi=[];end
if nargin<14;ssTi=[];end

if FLAGS0(1)==0;eval(['load ',fl_nr]);Atot=Apr;end
FLAGS=FLAGS0;SCALES=SCALES0;LSMN=LSMN0;cl_spec=cl_spec0;
hold off, format short e,
[nn,ninp]=size(u);[nn,noutp]=size(y);
valid=FLAGS(1);flag=FLAGS(2);
if length(FLAGS)>=3;ZER_MEAN=FLAGS(3);else;ZER_MEAN=0;end
stp=SCALES(1);tu_scal=SCALES(2:3);
fpts=SCALES(4);Mwin=SCALES(5);
nc=LSMN(1);red_meth=LSMN(2);
cutoff=LSMN(3:length(LSMN));
if length(cutoff)<5;cutoff(5)=1;end
if ZER_MEAN == 1
 %   VarU=sqrt(diag(var_nan(u))); VarY=sqrt(diag(var_nan(y)));
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
if ZER_MEAN >= 1
 %   VarU=sqrt(diag(var_nan(u))); VarY=sqrt(diag(var_nan(y)));
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
if ZER_MEAN <= -1
 %   VarU=sqrt(diag(var_nan(u))); VarY=sqrt(diag(var_nan(y)));
    nan_ind=[find(isnan(u(:,1)));nn+1];
    n_nan=length(nan_ind);f1_1=1;
    for i_nan = 1:n_nan
        ind_nan=[f1_1:nan_ind(i_nan)-1];
        O_V = ones(length(u(ind_nan,:)),1);
        u(ind_nan,:)=u(ind_nan,:)-O_V*u(f1_1,:);
        y(ind_nan,:)=y(ind_nan,:)-O_V*y(f1_1,:);
        f1_1=nan_ind(i_nan)+1;
    end
end

    if abs(ZER_MEAN) ==2 
            u=u*inv(VarU);
            y=y*inv(VarY);
    end 

if valid==0 
    scalinp=1;t_dil=1;
else
    scalinp=tu_scal(2);t_dil=tu_scal(1);stp=stp/t_dil;
    F_bw=F_bw*t_dil;cl_spec([1,3])=cl_spec([1,3])*t_dil;
    P_num=P_num*t_dil;P_den=P_den*t_dil;
    u=u*scalinp;
end
while length(F_ord)<noutp,
    F_ord=[F_ord,F_ord(length(F_ord))];
end
while length(F_bw)<noutp,
    F_bw=[F_bw,F_bw(length(F_bw))];
end
fil_BW=log10(min(F_bw));
frw=logspace(-2+fil_BW,2+fil_BW,100)';
if length(cutoff)<5;
    while length(cutoff)<5;cutoff=[cutoff,1];end
end
if (cutoff(1))>0;cutoff(1)=max(cutoff(1),1.e-4);end
if (cutoff(2))<1.e-4;cutoff(2)=1.e-4;end
if (cutoff(3))<1.e-4;cutoff(3)=1.e-4;end
if (cutoff(4))<1.e-4;cutoff(4)=1.e-4;end
if (cutoff(5))<1.e-4;cutoff(5)=1.e-4;end

%------------------------- 
t=[0:nn-1]'*stp;figure(1);clf;
subplot(121),plot(t,y(:,1));title('1st output')
subplot(122),plot(t,u(:,1));title('1st input')
uf=u;yf=y;
%--------q2y=1;del0=.25;norflag=0;r1=1;r2=1;

%-------------------------
if valid ~=0
    
    %---------------------------------------- filters 
    nuf=poly([-P_num]); def=poly([-P_den]);
    if P_num ~= P_den
        filternorm = norm(tf(nuf,def),'inf');nuf = nuf/filternorm;
        uf=mimofilt(nuf,def,u,t,IC_pts);yf=mimofilt(nuf,def,y,t,IC_pts);
    end
end
if m_red(1)==-1;uf1=uf;uf=yf;yf=uf1;disp('inverted model estimation');end

nan_ind=[find(isnan(uf(:,1)));length(uf)+1];
%--------------------------------------------- balancing
if valid ~=0
    CCC=[];FFF=[];num=[10];
    for i=1:noutp
        filt_i=ss(-F_bw(i),F_bw(i),1,0); for ii=2:F_ord(i),filt_i=filt_i*ss(-F_bw(i),F_bw(i),1,0);end
        n=F_ord(i);
  %      den=poly(-ones(1,F_ord(i))*F_bw(i));
  %      [f,q,cc,dd]=tf2ss(num,den);
        f=filt_i.a;q=filt_i.b;
        W=lyap(f,q*q');[U,S,V]=svd(W);TR=sqrt(inv(S))*U';
        f=TR*f*inv(TR);q=TR*q;
        CCC=[CCC,zeros(i-1,n);zeros(1,length(FFF)),q'];
        FFF=[FFF,zeros(length(FFF),n);zeros(n,length(FFF)),f'];
    end
    %----------------------------------------------------------
    disp('----------PARAMETER ESTIMATION----------------------')
    %----------------------------------------------------------
    if valid ==2
        %      par_estx
        [thx,err,nan_ind]=...
            par_estx(ninp,noutp,F_ord,FFF,CCC,uf,yf,t,cutoff,nc,red_meth);
        [Atot,Btot,Ctot,Dtot,Ftot,Bytot,Dytot,X0]=...
            mod_xtrx(FFF,CCC,thx,ninp,noutp,F_ord,nan_ind);
    elseif valid==3
        %      par_esta
        [thx,err,nan_ind]=...
            par_esta(ninp,noutp,F_ord,FFF,CCC,uf,yf,t,cutoff,nc,red_meth);
        [Atot,Btot,Ctot,Dtot,Ftot,Bytot,Dytot,X0]=...
            mod_xtra(FFF,CCC,thx,ninp,noutp,F_ord,nan_ind);
    else
        %      par_est
        [thx,err,nan_ind]=...
            par_est(ninp,noutp,F_ord,FFF,CCC,uf,yf,t,cutoff,nc,red_meth);
        [Atot,Btot,Ctot,Dtot,Ftot,Bytot,Dytot,X0]=...
            mod_xtr(FFF,CCC,thx,ninp,noutp,F_ord,nan_ind);
    end
 
    figure(2);clf;plot(t,err);title('Estimation Errors');
display(['estim. error std = ',num2str(sqrt(sum(var_nan(err))))  ])
display(['filt. input/output std = ',num2str(sqrt( sum(var_nan([uf  yf])))) ])
    pause(0.5)
%    bwestim(err,uf,yf);
 
    if flag>=0
        disp('Identified plant eigenvalues (scaled)')
        EAT=eig(Atot)
        disp('Identified plant transmission zeros (scaled)')
        TAT=tzero(Atot,Btot,Ctot,Dtot)
        if flag>1
            querytem=input('Gain Feedback Evals (eig(A-ans*B*C))? (y~=0) [0] ');
            if isempty(querytem);querytem=0;end
            if querytem ~= 0 & ninp==noutp
                FAT=eig(Atot-Btot*Ctot*querytem)
            end
        end
    end
    
    %    [frun,AUNC]=unc_bndf(err,u,Ftot,Bytot,Ctot,Dytot,3,fpoints,stp);
    %  [frun,AUNC]=unc_bn(err,u,u,1,1,fpoints,stp);
    %  ssdf=sv3_5(Ftot,Bytot,-Ctot,-Dytot,3,frun);
    %  AUNC=AUNC./(min(ssdf)');
    %  loglog(frun,AUNC);grid;title('Additive Uncertainty');pause
sspla=sv3_5(Atot,Btot,Ctot,Dtot,1,frw);
figure(3);clf;loglog(frw,sspla);title('Identified Plant Singular Values');
    
    %----------------------------------------------------------
    if m_red(1)~=0
        disp('----------MODEL REDUCTION----------------------')
        %----------------------------------------------------------
        [Apr,Bpr,Cpr,Dpr,Fpr,Bypr,Dypr]=...
            mod_red(Atot,Btot,Ctot,Dtot,Ftot,Bytot,Dytot,F_bw,flag,m_red);
    else
        Apr=Atot;Bpr=Btot;Cpr=Ctot;Dpr=Dtot;Fpr=Ftot;Bypr=Bytot;Dypr=Dytot;
        MR=[];
    end
Bpr=Bpr*VarU; Cpr=VarY*Cpr; Dpr=VarY*Dpr*VarU;

end
%----------------------------------
u=u/scalinp;uf=uf/scalinp;Bpr=Bpr*scalinp;Dpr=Dpr*scalinp;%MR=MR*scalinp

sspla=sv3_5(Apr,Bpr,Cpr,Dpr,1,frw);
figure(3);clf;loglog(frw,sspla);title('Identified Plant Singular Values');
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
if size(X0)==size(X0i)
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
    if flag>1
        querytem=input('Gain Feedback Evals ? (y~=0) [0]  ');
        if isempty(querytem);querytem=0;end
        if querytem ~= 0 & ninp==noutp
            FAT=eig(Apr-Bpr*Cpr*querytem)
        end
    end
end

disp('Identified DC gain')
DC=Dpr-Cpr*inv(Apr)*Bpr
%----------------------------------------------------------
disp('----------IDENTIFICATION TESTS----------------------')
%----------------------------------------------------------

%  querytem=input('Perform frequency domain analysis ? (1=y) [0] ');
%     if querytem ==1
%        FREQan
%     end

if flag>=1
    querytem=input('Compare output time responses? (0=no) [1]  ');
    if isempty(querytem);querytem=1;end
else
    querytem=0;
end
Y=lsim_nan(Apr,Bpr,Cpr,Dpr,u,t,X0);
display(['Add. error std = ',num2str(sqrt(sum(var_nan(Y-y))))  ])
display(['input/output std = ',num2str(sqrt( sum(var_nan([u  y])))) ])
if querytem ~=0
    for i=1:noutp
        figure(4);clf
        plot(t,y(:,i),'r',t,Y(:,i),'g');
        title('Actual (r) vs. Identified (g) system comparison');
        pause
    end
end
figure(4);clf
plot(t,y-Y);title('Additive Identification Error');
pause(0.5)

%----------------------------------------------------------
disp('----------UNCERTAINTY ESTIMATION----------------------')
%----------------------------------------------------------
frun=[];AUNC=[];MUNC=[];succ=1;
if flag>=1
    querytem=input('Perform uncertainty analysis (UNCERF)? (0=no) [1]  ');
    if isempty(querytem);querytem=1;end
else
    querytem=1;
end
if querytem ~=0
    % ----------------     UNCERF
    [frun,AUNC,MUNC,succ,cl_spec]=uncerf(Apr,Fpr,Bpr,Cpr,Dpr,Bypr,Dypr,X0,u,y,t,fpts,flag,cl_spec,ssTi,ssDmi);
end
%------------------------------- recover time dilation
Apr=Apr/t_dil;Bpr=Bpr/t_dil;Fpr=Fpr/t_dil;Bypr=Bypr/t_dil;stp=stp*t_dil;
frw=frw/t_dil;frun=frun/t_dil;t=t*t_dil;DC=Dpr-Cpr*inv(Apr)*Bpr;
F_bw=F_bw/t_dil;cl_spec([1,3])=cl_spec([1,3])/t_dil;
P_num=P_num/t_dil;P_den=P_den/t_dil;

%------------------------------- save results

FLAGS=[valid,flag,ZER_MEAN];tu_scal=[t_dil,scalinp];
SCALES=[stp,tu_scal,fpts,Mwin];
LSMN=[nc,red_meth,cutoff];
SYS=ss(Apr,Bpr,Cpr,Dpr);

figure(8)
step(SYS);



if flag >=2 | length(fl_ns)<2
    savsys=input(['Give filename to save ID results; current: ',fl_ns,'  '],'s');
else
    savsys=fl_ns;
end

if length(savsys)>=2
    eval(['save ',savsys,' Apr Bpr Cpr Dpr Fpr Bypr Dypr X0 DC frw sspla MR frun MUNC AUNC F_bw F_ord cl_spec P_num P_den IC_pts m_red LSMN SCALES FLAGS'])
    f_sav=savsys;
else
    f_sav=[];
end
