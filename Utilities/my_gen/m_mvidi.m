function succ=m_mvidi(ti,ui,yi,fl_nr,fl_ns,flag);
% function succ=
%    m_mvidi(ti,ui,yi,fl_nr,fl_ns,flag);
%
%   Shell/calling function for system identification
%   ti (opt.), ui, yi are the I/O pairs at ti
%   fl_nr= optional read file to establish defaults (from previous ID)
%   fl_ns= optional write file
%   bw = intended closed-loop bandwidth or auxiliary filter BW.
%   f_ord = filter order
%   rolloff = intended rolloffs for cl.lp. T,S sensitivities
%   flag= 0 for auto, 1 or higher for intermediate inputs and displays
%   
%     Comments:  flag=0 should normally be used after a preliminary ID
%   or when the ID bandwidth has been previously established.
%                Negative values for bw set the filter bandwidth; if a 
%   vector, then each component becomes the channel bandwidth. 
%   If bw is positive, it is interpreted as intended closed-loop BW
%   and all channels are ID'd with the same BW filter; 
%   in this case, it can be a vector and its second element is taken 
%   as the sensitivity break frequency.
%                f_ord can be a vector specifying the order of the
%   auxiliary ID filters per channel.
%                rolloff specifies the rolloff rates for T and S in
%   the same order (default [2,2], i.e., "integrating" plants,
%   or [rolloff,1] if scalar).  rolloff can be a 4-element vector
%   with the last two entries specifying RISK and AGRS parameters.
%                The rolloff and the closed-loop BW are used in the
%   uncertainty computations and, eventually, as closed-loop specs; 
%   Notice that a rough estimate can be used that can be altered 
%   during the controller design phase since the uncertainty 
%   computations are relatively incensitive to the exact values. 
%   The function "unc_chk" can be used after the controller design 
%   to establish the confidence of the design in a more optimal 
%   manner (needed only for tight specs).
%                More parameters (prefilters, LS/MN, etc.) can be 
%   specified interactively via the menu screen.  The step size 
%   is automatically computed from ti but can be overridden later.  
%   Also, if a time dilation is specified, (e.g., to obtain models 
%   with more favorable pole-zero structure) all specs should be 
%   given in the original time scale. Although the intermediate ID 
%   displays are given in the dilated time scale, the final results
%   are converted and saved in the original time scale.
%
%   ssTi,ssDmi = inner loop complementary sv's and uncertainty (rows)
%

% KST 8/97   

%--- Default Settings -----------------------
bw_fact=1;bw_cl=5;TS_roll=[2 2];F_ord=3;F_bw=0;
if length(ti)>2;stp=ti(2)-ti(1);elseif length(ti)==1;stp=ti;else;stp=1/60;end
P_num=1;P_den=1;IC_pts=1;tu_scal=[1 1];valid=1;
nc=0;cutoff=[.2 2 2 2 1];red_meth=2;
fpts=80;Mwin=4;t_dil=1;m_red=[0 3 0];
RISK=1;AGRS=1;
FLAGS=[];SCALES=[];LSMN=[];ZER_MEAN=0;
[ntot,ninp]=size(ui);[ntot,noutp]=size(yi);
U_SCAL=ones(ninp,1);Y_SCAL=ones(noutp,1);
targ_T=[];targ_S=[];
%--------------------------------------------
if nargin <4;fl_nr=[];end
if nargin <5;fl_ns=[];end
if nargin<6;flag=1;end;if isempty(flag);flag=1;end
%------------------load file 
          if length(fl_nr)>=2;
             eval(['load ',fl_nr]);
               if length(SCALES)>0;  %----------------new format
                  valid=FLAGS(1);
                  if length(FLAGS)>=3;ZER_MEAN=FLAGS(3);else;ZER_MEAN=0;end
                  stp=SCALES(1);tu_scal=SCALES(2:3);fpts=SCALES(4);
                  Mwin=SCALES(5);
                  nc=LSMN(1);red_meth=LSMN(2);cutoff=LSMN(3:length(LSMN));
                  if length(cutoff)<5;cutoff(5)=1;end
               else                  %----------------old format
                  tu_scal(1)=t_dil;F_bw=mean(abs(roots(den)));
                  bw_cl=F_bw*bw_fact;F_ord=length(den)-1;
               end
          end
%---------------------------

%--- Interactive specs -----------------------
%---------------------------------------------------
if flag >0
z_fig=figure;
   if z_fig>1
     close(z_fig)
   else
disp('**********************************************************')
disp('!!!!!!!   PLEASE ADJUST FIGURE WINDOW SIZE    !!!!!!!!!!!!')
disp('**********************************************************')
pause
   end
end
%----------------------------------------------------
if ~(flag<=0 & length(fl_nr)>1)
   temp=1;
   while temp ~=0
     pause(0.5);home
     disp('      INPUT OPERATION MENU')
     disp('------ 1. Closed-Loop Specs')
     disp('------ 2. Auxiliary Filters (Order and Bandwidth)')
     disp('------ 3. Sampling Time, Time dilation and input scaling')
     disp('------ 4. Fitting parameters: weighted LS')
     disp('------ 5. Prefiltering')
     disp('------ 6. Estimation Constraints')
     disp('------ 7. Model order reduction')
     disp('------ 8. Trace flag, estimation/validation flag')
     disp('------ 9. Read/write data from/to a file')
     disp('------ 0. Exit Menu and Begin')
     disp(' ')
     temp=input('---Select Input Operation      ');
     if length(temp)~=1;temp=1;end
%       ----------------------------------------  
        if temp==1
          str1=[num2str(bw_cl)];
          disp(['Current defaults: cl.lp. bw     [',str1,']'])
          temp2=input('Enter desired closed-loop bandwidth   ');
            if length(temp2)== 1;
              bw_cl=temp2;
            end
%       ----------------------------------------  
        elseif temp==2
          disp('Current default values: filter orders and bandwidth')
          disp(F_ord),disp(F_bw)
          temp2=input('Enter aux. filter orders   ');
          if length(temp2)>=1;F_ord=temp2;end
          temp2=input('Enter aux. filter bandwidth (0=auto from cl_bw)   ');
          if length(temp2)>=1;F_bw=temp2;end
          if F_bw==0;F_bw=bw_cl/bw_fact;end
%       ----------------------------------------  
        elseif temp==3
          disp('Current default values: sampling time and t/u scaling')
          disp([stp,tu_scal])
          temp2=input('Enter sampling time   ');
          if length(temp2)<1;temp2=1;else;stp=temp2;end
          temp2=input('Enter [time dilation, input scaling]   ');
          if length(temp2)<1;temp2=1;
            elseif length(temp2)==1;tu_scal(1)=temp2;
            elseif length(temp2)==2;tu_scal=temp2;
          end
%       ----------------------------------------  
        elseif temp==4
          disp(['LS/MN estimation defaults: cutoff, stability, diagonal, feedback, IC weights  ']),disp(cutoff)
          temp2=input('cutoff parameters  ');
          if length(temp2)>=1;
             cutoff(1:length(temp2))=temp2;
          end
          temp2=input(['method (1=svd, 2=MN/LS); current:  ',num2str(red_meth),'  ']);
          if length(temp2)==1;red_meth=temp2;end
%       ----------------------------------------  
        elseif temp==5
          disp('Prefilter numerator and denominator polynomials')
          disp(P_num),disp(P_den) 
          temp2=input('numerator (-zeros)   ');
          if length(temp2)>=1;P_num=temp2;end
          temp2=input('denominator (-poles)  ');
          if length(temp2)>=1;P_den=temp2;end
          temp2=input(['Initial condition averaging points  [', num2str(IC_pts),']  ']);
          if length(temp2)>=1;IC_pts=temp2;end
          disp('Input/Output scaling factors')
          disp(U_SCAL');disp(Y_SCAL')
          temp2=input('Input scales  ');
          if length(temp2)>=1;U_SCAL(1:length(temp2))=temp2';end
          temp2=input('Output scales  ');
          if length(temp2)>=1;Y_SCAL(1:length(temp2))=temp2';end

%       ----------------------------------------  
        elseif temp==6
          temp2=input(['Constraints: 1=strictly proper, 0=biproper;  current:  ',num2str(nc),'  ']);
          if length(temp2)>=1;nc=temp2;end
%       ----------------------------------------  
        elseif temp==7
          str1=[num2str(m_red(1)),',',num2str(m_red(2)),',',num2str(m_red(3))];

          disp(['Model Order Reductions defaults:  [',str1,']'])
          disp('1st=enable switch, 2nd=method(1=Schur, 2=relative, 3=interactive)')
          disp('3rd = threshold (Schur: absolute, Relative: dB)') 
          temp2=input('Model order reduction parameters   ');
          if length(temp2)>=1;m_red(1:length(temp2))=temp2;end
%       ----------------------------------------  
        elseif temp==8
          disp('Trace/Input flag: higher values=more intermediate outputs')
          temp2=input(['Trace/input flag: current  ',num2str(flag),'   ']);
          if length(temp2)>=1;flag=temp2;end
          disp('Estimation (=1,2) /Validation (=0) flag')
          temp2=input(['Estimation/Validation flag: current  ',num2str(valid),'   ']);
          if length(temp2)>=1;valid=temp2;end
          temp2=input(['Extract mean from data? (1=y) [', num2str(ZER_MEAN),']  ']);
          if length(temp2)>=1;ZER_MEAN=temp2;end
%       ----------------------------------------  
        elseif temp==9
          temp2=input(['Read data from file;  current:  ',fl_nr,'   '],'s');
          if length(temp2)>=2;
             fl_nr=temp2; eval(['load ',fl_nr]);
               if length(SCALES)>0;  %----------------new format
                  valid=FLAGS(1);
                  if length(FLAGS)>=3;ZER_MEAN=FLAGS(3);else;ZER_MEAN=0;end
                  stp=SCALES(1);tu_scal=SCALES(2:3);fpts=SCALES(4);
                  Mwin=SCALES(5);
                  nc=LSMN(1);red_meth=LSMN(2);cutoff=LSMN(3:length(LSMN));
			if length(cutoff)<5;cutoff(5)=1;end
               else                  %----------------old format
                  tu_scal(1)=t_dil;F_bw=mean(abs(roots(den)));
                  bw_cl=F_bw*bw_fact;F_ord=length(den)-1;
               end
          end
          temp2=input(['Save data to file;  previous:  ',fl_ns,'   '],'s');
          fl_ns=temp2;
%       ----------------------------------------  
        elseif temp==0
          disp('Good Luck!')
%       ----------------------------------------  
        else
          why
          pause(2)
%       ----------------------------------------  
        end
   end
end

%--- Call function ---------------------------
FLAGS=[valid,flag,ZER_MEAN];
SCALES=[stp,tu_scal,fpts,Mwin];
LSMN=[nc,red_meth,cutoff];
          if F_bw==0;F_bw=bw_cl(1)/bw_fact;end

[xx,succ]=m_mvid(ui,yi,F_bw,F_ord,bw_cl,P_num,P_den,IC_pts,m_red,fl_nr,fl_ns,LSMN,SCALES,FLAGS,U_SCAL,Y_SCAL,targ_T,targ_S);

if length(xx)>=2 
   disp(['***  ID data and results saved in file  ',xx])
else
   disp(['***  ID data not saved !!'])
end
if max(succ)<1
   disp('***  ID seems successful for the desired closed-loop specs.')
   semilogy(succ);title('Closed-loop Robust Stability test');
   xlabel('SUCCESS (?)');pause(1)
else
   disp('*** ID may be UNSUCCESSFUL for the desired closed-loop specs !!')
   semilogy(succ);title('Closed-loop Robust Stability test');
   xlabel('FAILURE (?)');pause(1)
end

