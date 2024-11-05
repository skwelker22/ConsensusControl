function [A,B,C,D]=mh_wgt(nout,tar);
% [A,B,C,D]=mh_wgt(nout,tar);
% For the interactive specification of frequency weights (Hinf)
% The weight is a diagonal of bi-proper transfer functions
%     defined in terms of poles, zeros &  DC gain in tar
% Upon completion the inverse t.f.'s are realized in s.s., 
%    concatenated and returned in [A,B,C,D]

if nout~=tar(1);disp('mh_wgt: incorrect dimensions');return;end

ind=1;
for ii=1:nout
    den=1;num=1;nord=tar(ind+1);
    w0=tar(ind+2:ind+1+nord);w1=tar(ind+2+nord:ind+2*nord+1);r1=tar(ind+2+2*nord);
    for i=1:nord;den=conv(den,[1/w0(i),1]);end
    for i=1:nord;num=conv(num,[1/w1(i),1]);end
    num=real(num);den=real(den);
    [A1,B1,C1,D1]=tf2ss(den,r1*num);
      if ii==1
        A=A1;B=B1;C=C1;D=D1;
      else
        [A,B,C,D]=append(A,B,C,D,A1,B1,C1,D1);
      end
    ind=ind+2+2*nord;
end



