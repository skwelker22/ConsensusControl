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
n1=1;n2=nan_ind(1)-1;
s_e=0*frun;s_u=0*frun;s_y=0*frun;
for i=1:n_nan
   ui=u(n1:n2,:);yi=y(n1:n2,:);ei=err(n1:n2,:);
  % ui=ui.*hanning(length(ui));yi=yi.*hanning(length(ui));ei=ei.*hanning(length(ui));
   nn=n2-n1+1;nn2=round(nn/2);
   frunx=fft_freq(nn,stp);frunx=frunx(1:nn2);
   fer=fft(ei);
   for ii=1:noutp,fer(:,ii)=fer(:,ii)-ei(1,ii)/2-ei(nn,ii)/2;end
   fer=abs(fer(1:nn2,:));
   fer=(fer.*fer)*ones(noutp,1);fer=sqrt(fer);
   fu=fft(ui);
   for ii=1:ninp,fu(:,ii)=fu(:,ii)-ui(1,ii)/2-ui(nn,ii)/2;end
   fu=abs(fu(1:nn2,:));
   fu=(fu.*fu)*ones(ninp,1);fu=sqrt(fu);
   fy=fft(yi);
   for ii=1:nouty,fy(:,ii)=fy(:,ii)-yi(1,ii)/2-yi(nn,ii)/2;end
   fy=abs(fy(1:nn2,:));
   fy=(fy.*fy)*ones(nouty,1);fy=sqrt(fy);

   f_ind=find(fu > toler & fy > toler);
   frunx=frunx(f_ind);fer=fer(f_ind);fu=fu(f_ind);fy=fy(f_ind);
   fer=smoothin(frunx,fer,frun);
   fu=smoothin(frunx,fu,frun);
   fy=smoothin(frunx,fy,frun);
   s_e=s_e+fer;s_u=s_u+fu;s_y=s_y+fy;
end

