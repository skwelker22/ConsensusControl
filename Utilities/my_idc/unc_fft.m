function [frun,s_e,s_u,s_y]=unc_fft(err,u,y,pts,stp);
%function [frun,s_e,s_u,s_y]=unc_fft(err,u,y,pts,stp);
% computes the fft bounds for uncertainty estimation
% in sys.id. (used by MIMOID/UNCERF)

toler=1.e-8;
[nnt,noutp]=size(err);[nnt,ninp]=size(u);
[nnt,nouty]=size(y);
i_nnan=find(~isnan(u(:,1)));
nan_ind=[find(isnan(u(:,1)));nnt+1];
nan_ind1=[1;find(isnan(u(:,1)))];
nn_s=max(nan_ind-nan_ind1);
n_nan=length(nan_ind);

if length(pts)==1
  frmin=2*pi/nn_s/stp;frmax=pi/stp/4;
  frun=logspace(log10(frmin),log10(frmax),pts)';
  frun2=logspace(log10(frmin),log10(frmax),2*pts)';
else
  frun=pts;
end

%  ---------------------------------------- perform ffts
n1=1; s_e=0*frun;s_u=0*frun;s_y=0*frun;
for i_nan=1:n_nan
   ind_i=[n1:nan_ind(i_nan)-1]; n1=nan_ind(i_nan)+1;
   ui=u(ind_i,:);yi=y(ind_i,:);ei=err(ind_i,:);
  % ui=ui.*hanning(length(ui));yi=yi.*hanning(length(ui));ei=ei.*hanning(length(ui));
  % pwelch(..,2pi/stp)*N/2pi = pwelch(.)*N*stp/2pi = (fft*stp)^2 = F()^2
   nn=length(ind_i);nn2=round(nn/2);

   [EI,frw]=pwelch(ei(:,1),[],[],[],2*pi/stp);
   scf=length(frw)*stp/2/pi;
   EI=pwelch(ei(:,1))*scf;UI=pwelch(ui(:,1))*scf;YI=pwelch(yi(:,1))*scf;
   for i=2:noutp;EI=EI+pwelch(ei(:,i))*scf;end
   for i=2:nouty;YI=YI+pwelch(yi(:,i))*scf;end
   for i=2:ninp;UI=UI+pwelch(ui(:,i))*scf;end
   
   fer=smoothin(frw,max(toler,EI),frun);
   fu=smoothin(frw,max(toler,UI),frun);
   fy=smoothin(frw,max(toler,YI),frun);
   
   s_e=s_e+fer;s_u=s_u+fu;s_y=s_y+fy;
end

s_e=sqrt(s_e);s_u=sqrt(s_u);s_y=sqrt(s_y);