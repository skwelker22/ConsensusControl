function [h,phi]=PIDF_OBJ(W,zi,K,op_t,fr);
% function [h,phi]=PIDF_OBJ(W,zi,K,op_t,fr);
% Defines the optimization objective in pidftune

% K. TSAKALIS, 8/23/96


temp = W*K-zi;npt=length(temp);

if op_t(1)==1
    dw   = [fr;0]-[0;fr];dw=dw(1:npt);
    ph   = (abs(temp)).^2;
    phi  = (ph(1:npt-1)+ph(2:npt))'*dw(2:npt)/2+ph(1)*dw(1);
    h    = (temp(1:npt-1).*dw(2:npt))'*W(1:npt-1,:);
    h    = h+(temp(2:npt).*dw(2:npt))'*W(2:npt,:);
    h    = real(h+2*temp(1)'*W(1,:)*dw(1));
    h    = h'/norm(h);
elseif op_t(1) < 0
    phim = max(abs(temp))^2;
    phi  = max(abs(temp))^2 - op_t(1)*op_t(1);
    if phi > 0
        indp = min(find(abs(temp)==sqrt(phim)));
        h    = 2*(real(temp(indp)'*W(indp,:)))';
    else
        phi = 0; h = 0*K;
    end
else
    phi  = max(abs(temp))^2;
    indp = min(find(abs(temp)==sqrt(phi)));
    h    = 2*real(temp(indp)'*W(indp,:));
    h    = h'/norm(h);
end
