% Script file FREQAN, for the frequency domain analysis
%   of system id from data; used by MIMOID
%   Requires I/O pairs [u,y] to exist in the workspace
%      stp (sampling time), ninp,noutp (number of inputs/outputs)
%   For comparisons with parametric id, the identified frequency
%      response must also be in the workspace (frw,MR)
%   Convention:  MR... contains freq.resp. in colums 
%   Options: generic FFT id and Ljung's spectral analysis
%   Outputs: [fftfr,Gest] (FFT), [frg,Gspa,SDGspa] (Ljung)

querytem=input('Perform frequency domain estimation ? (1=y) ');
   if querytem ==1
     fpts=8;nover=4;M=8;
     if nn<=2^11;fpts=4;M=4;end
     if nn<=2^9;fpts=2;M=2;end
     if nn<=2^7;fpts=1;M=2;nover=2;end
     win_pts=fix(nn/M);over_pts=fix(win_pts/nover);
     fr_pts=fix(win_pts/fpts);

     Gest=zeros(length([1:fr_pts]),ninp*noutp);
       for i=1:noutp
         [fftfr,Gest(:,i:noutp:ninp*noutp)]=...
           FREQmdl(y(:,i),u,stp,M,nover,fpts);
       end
   clg
     for i=1:ninp
       ind_set=[(i-1)*noutp+1:i*noutp];
       loglog(frw,abs(MR(:,ind_set)),'-',...
             fftfr,abs(Gest(:,ind_set)),'--');
       pause
     end
   end

querytem=input('Perform spectral analysis ? (1=y) ');
   if querytem ==1
   frg=[2/nn:(1/4-2/nn)/255:1/4]'*pi/stp;
   Gspa=zeros(length(frg),ninp*noutp);SDGspa=Gspa;
       for i=1:noutp
         [Gpla,NS]=spa([y(:,i), u],nn/8,frg',-1,stp);
         ind_set=find(Gpla(1,:)<10)
         Gspa(:,i:noutp:ninp*noutp)=Gpla(2:length(Gpla),ind_set);
      SDGspa(:,i:noutp:ninp*noutp)=Gpla(2:length(Gpla),ind_set+1);
       end
clg

     for i=1:ninp
       ind_set=[(i-1)*noutp+1:i*noutp];
       loglog(frw,abs(MR(:,ind_set)),'-',...
         frg,abs(Gspa(:,ind_set)),'--');
       pause
     end
   end
