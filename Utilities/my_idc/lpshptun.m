function [lshape,lam,al,pm]=lpshptun(frw,lshapeo,lam,al,pm);
% function [lshape,lam,al,pm]=lpshptun(frw,lshapeo,lam,al,pm);
%   Interactive specification of a 2nd-order loop-shape 
%               lam*[s al]/[s^2 + pm s] 
%   for PID tuning

clg;
rdone=0;lshape=lshapeo;
while rdone==0
  disp('Target closed loop poles')
  roots([1,pm+lam,lam*al])
  loglog(frw,(1)./abs(1+lshape),frw,abs(lshape./(1+lshape)))
  title('Target sensitivities')
  pause
  T=[0:.1/lam:10/lam]';
  Y=step(lam*[1 al],[1 lam+pm lam*al],T);
  plot(T,Y);title('Target Step Response')
  pause
  rdone=input('Done?  ');
    if rdone == 0
      mat=input('Give adjustment multipliers for lam,al,pm  ')
      lam=lam*mat(1);al=al*mat(2);pm=pm*mat(3);
      lshape=freqs(lam*[1 al],[1 pm 0],frw).'; 
    end 

end
