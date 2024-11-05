function [A,B,C,D]=hinf_wgt(n1,w0,w1,r1);
% [A,B,C,D]=hinf_wgt(n1,w0,w1,r1);
% For the interactive specification of frequency weights (Hinf)
% The weight is a diagonal of a bi-proper transfer function
%    (n1 replicas) defined in terms of poles, zeros &  DC gain 
% Upon completion the inverse t.f.'s are realized in s.s., 
%    concatenated and returned in [A,B,C,D]

if r1~=0
nw1=0;den=1;num=1;
    for i=1:length(w0);den=conv(den,[1/w0(i),1]);end
    for i=1:length(w1);num=conv(num,[1/w1(i),1]);end
    num=real(num);den=real(den);
[A1,B1,C1,D1]=tf2ss(den,r1*num);nw1=length(A1);
else
[A1,B1,C1,D1]=tf2ss(w0,w1);nw1=length(A1);
end

A=zeros(n1*nw1,n1*nw1);
B=zeros(n1*nw1,n1);
C=zeros(n1,n1*nw1);
D=zeros(n1,n1);
  for i=1:n1
    A(nw1*(i-1)+1:nw1*i,nw1*(i-1)+1:nw1*i)=A1;
    B(nw1*(i-1)+1:nw1*i,i)=B1;
    C(i,nw1*(i-1)+1:nw1*i)=C1;
    D(i,i)=D1;
  end
