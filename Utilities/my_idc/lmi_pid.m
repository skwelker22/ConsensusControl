function [K,Ub]=LMI_PID(W,zi,toler,T,pos,B,BDZ);
% function [K,Ub]=LMI_PID(W,zi,toler,T,pos,B,BDZ);
%  Convex optimization function to perform PID tuning
%  Calls to LMI_upd and PID_cns
%
K=[0;0;0]; Aell=B(1)*B(1)*eye(length(K),length(K));
L=[];U=[];err=1;erro=1;kount=0;phi=0;kkount=0;indp=0;

while err>toler
  [h,psi]=PID_cns(K,1.e-4,T,pos,B(2),BDZ);h0=norm(h);
    if h0>0
      % 'constraint-iteration'
      iterty=-1;deepcut=psi;
    else
      % 'objective-iteration'
      iterty=1; temp=W*K-zi;
      phi=max(abs(temp))^2;indp=min(find(abs(temp)==sqrt(phi)));
      h=2*real(temp(indp)'*W(indp,:));h=h'/norm(h);deepcut=0;
      erro=sqrt(h'*Aell*h);L=[L;phi-erro];U=[U;phi];
    end
[K,Aell]=LMI_upd(K,h,Aell,deepcut);
kount=kount+1;kkount=kkount+1;
h0=norm(PID_cns(K,1.e-4,T,pos,B(2),BDZ));err=(h0)+erro;
    if kkount>=20,kkount=0;
%      disp([kount, err, phi,iterty,indp])
    end
end
Ub=sqrt(min(U));


