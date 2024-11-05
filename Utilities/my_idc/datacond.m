function [ys,TS]=datacond(y,s_ty);
% function [ys,TS]=datacond(y,s_ty);
%   data conditioning of column vector data y using covariance matrix transformation
%   s_ty: 3-d row vector defining the transformation attributes:
%         s_ty(1)=1: extract mean, 0: not
%         s-ty(2)=0: full covariance inverse, 1: diagonal only
%         s_ty(3)=1: recover mean, 0: not, -1: match 1st row

if nargin < 2; s_ty=1;end
while length(s_ty)<3;s_ty=[s_ty,0];end
Ix=ones(length(y),1);
YM=mean(y);
Y=y-Ix*YM*s_ty(1);
YC=Y'*Y/length(y);
TS=eye(size(YC));
[U,S,V]=svd(YC);
if min(diag(S))<1e-15
   disp('datacond: covariance matrix near singular')
   diag(S)
else
   TS=inv(sqrt(S))*U';
end
if s_ty(2)==1
  TS=diag(diag(TS));
end
ys=Y*TS';yms=mean(ys);
if s_ty(3)==1
  ys=ys+Ix*(YM-yms);
elseif s_ty(3)<0
  ys=ys+Ix*(y(1,:)-ys(1,:));
end
