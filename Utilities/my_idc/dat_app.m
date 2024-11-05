function [tx,Upx,Ysx,Ypx,Ywx]=dat_app(nx,zerom,scales);
%function [tx,Upx,Ysx,Ypx,Ywx]=dat_app(nx,zerom,scales);
%appends data files

done=0;first=0;
scal=ones(1,nx);
if length(zerom)==1;zerom=zerom*scal;end
if nargin>=3
  if ~isempty(scales);
     scal(1:min(length(scales),nx))=scales(1:min(length(scales),nx));
  end
end

for i=1:nx
  Ys=[];Yp=[];Yw=[];Up=[];
  fnm = input(' filename  ','s');
  eval(['load ',fnm])
  ov=0*t+1;
  if zerom(i)>1
    Up=Up-ov*mean(Up(1:zerom(i),:));
    Yp=Yp-ov*mean(Yp(1:zerom(i),:));
    Ys=Ys-ov*mean(Ys(1:zerom(i),:));
  elseif zerom(i)==1
    Up=Up-ov*(Up(1,:));
    Yp=Yp-ov*(Yp(1,:));
    Ys=Ys-ov*(Ys(1,:));
  end
    Up=Up*scal(i);
    Yp=Yp*scal(i);
    Ys=Ys*scal(i);

    if first==0
      tx = t; Upx = Up; Ysx = Ys; Ypx = Yp; Ywx = Yw;
    else
      tx = [tx;nan;t+first];
      Upx = [Upx;Up(1,:)*nan;Up];
      Ysx = [Ysx;Ys(1,:)*nan;Ys];
      Ypx = [Ypx;Ys(1,:)*nan;Yp];
      Ywx = [Ywx;Ys(1,:)*nan;Yw];
    end
  first = max(t)+1+first;
end
