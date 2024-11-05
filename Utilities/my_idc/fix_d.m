function D=fix_d(D,tol);
% function D=fix_d(D,tol);
% ensure that the min singular value of a matrix D
% is at least tol.

[UD,SD,VD]=svd(D);  [rd,cd]=size(D); nsd1=min(rd,cd);
SD1=diag(max(diag(SD),tol));
SD(1:nsd1,1:nsd1)=SD1(1:nsd1,1:nsd1);
D=UD*SD*VD';

