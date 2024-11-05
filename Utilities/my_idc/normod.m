function [F,C,Bu,Du,By,Dy]=normod(At,Bt,Ct,Dt,shift);

%function [F,C,Bu,Du,By,Dy]=normod(At,Bt,Ct,Dt);
% nonminimal realization of [At,Bt,Ct,Dt] based on a
% normalized coprime factorization (A,B,C,D);

[noutp,ninp]=size(Dt);
if nargin<5;shift=0;end

[AN,BN,CN,DN,AD,BD,CD,DD]=norm_cop(At,Bt,Ct,Dt,shift);

F=AN;Bu=BN;By=-BD;C=CN;Du=DN;Dy=eye(noutp,noutp)-DD;
