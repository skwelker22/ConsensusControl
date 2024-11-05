function [m,num,den]=awgt_sel(frw,clbw,clor,s_ty);
% [m,num,den]=awgt_sel(frw,clbw,clor,s_ty);
% frequency weight selection for ID/H-infty design.
% frw=frequency pts, clbw=closed-loop (S/T) BW, clor= (S/T) rolloff
% s_ty=0 for T.  DC/HF gains are fixed


if nargin < 4, s_ty=0;end

w0=ones(1,clor)*clbw;
w1=0*w0+1;
if s_ty==0;r1=1.4;else;r1=-1.7;end
if r1>0;w1=w1*max(w0)*300;else;w1=w1*min(w0)/200;end
den=1;for i=1:length(w0);den=conv(den,[1/w0(i),1]);end
num=1;for i=1:length(w1);num=conv(num,[1/w1(i),1]);end
if r1<0, r1=abs(r1)*den(1)/num(1);end
num=r1*num;
m=freqs(num,den,frw);
m=abs(m);
