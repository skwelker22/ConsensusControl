function [h,phi]=PIDqOBJ(W,zi,K);
% function [h,phi]=PIDqOBJ(W,zi,K);
% Defines the optimization objective in pidftune

% K. TSAKALIS, 8/11/04


      temp = W*K-zi;npt=length(temp);

      phi  = max(abs(temp))^2;
      indp = min(find(abs(temp)==sqrt(phi)));
      h    = 2*real(temp(indp)'*W(indp,:));
      h    = h'/norm(h);
