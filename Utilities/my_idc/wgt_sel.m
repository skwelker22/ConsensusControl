function [w0,w1,r1,frw,m]=wgt_sel(n1,fr,UN,s_ty,def_flag);
% [p,z,g,f,m]=wgt_sel(n1,fr,UN);
% Interactive frequency weight selection for H-infty design.
% The transfer function is defined in terms of poles, zeros &
%    DC or HF gain and the inverse magnitude is plotted 
%    against 1/UN (if specified).

noUN=1; if nargin<2;fr=0;end
if nargin <2 | length(fr)<1 | fr ==0
  frw=logspace(-3,3,100);noUN=0;else frw=fr;
end
if nargin < 3 | length(UN)~= length(frw), noUN=0;end
if nargin < 4, s_ty=0;end
if nargin < 5, def_flag=0;end


itemp=0;nw1=0;nw2=0;w1x=[];
while itemp==0
  nord=input('Weight order [2]   ');
  if length(nord)<1;nord=2;end
  w0=input('pole corner frequencies   ');
  if length(w0)<1;w0=1;end
    while length(w0)<nord;w0=[w0,w0(length(w0))];end
    if def_flag>=1
      w1x=input('zero corner frequencies   ');
      if length(w1x)<1;w1=1;else;w1=w1x;end
      while length(w1)<nord;w1=[w1,w1(length(w1))];end
    else
      w1=0*w0+1;
    end
    if def_flag>=1
      r1=input('DC/(-HF) gain (e.g., +1.4/-1.7)  ');
        if length(r1)<1;
         if s_ty==0;r1=1.4;else;r1=-1.7;end
        end
    else
      if s_ty==0;r1=1.4;else;r1=-1.7;end
    end
    if def_flag <1 | length(w1x)<1
      if r1>0;w1=w1*max(w0)*300;else;w1=w1*min(w0)/500;end
    end
    den=1;for i=1:length(w0);den=conv(den,[1/w0(i),1]);end
    num=1;for i=1:length(w1);num=conv(num,[1/w1(i),1]);end
    if r1<0, r1=abs(r1)*den(1)/num(1);end
  m=freqs(r1*num,den,frw);
    if noUN==0
      loglog(frw,abs(m));grid;
    else
      loglog(frw,sat(abs(m),40,1.e-4),...
            fr,sat((1)./abs(UN),40,1.e-4),'--');grid;
    end
  title('Target Sensitivities and uncertainty plots');pause
  itemp=input('done? (1=yes) [1]  ');
  if isempty(itemp);itemp=1;end
end

m=abs(m);
