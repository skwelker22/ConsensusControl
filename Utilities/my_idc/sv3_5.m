function [sv] = sigma(a,b,c,d,Type,w)
%
% Singular Value Frequency Response.
%
% SV = SIGMA(A, B, C, D, TYPE, W) produces the matrix SV
% containing the singular values of the square system : 
%                .
%                x = Ax + Bu
%                y = Cx + Du
%                                                   -1
%     with the frequency response G(jw) = C(jwI - A)  B + D
%
% SIGMA calculates the SVD of one of the following types:
%
%     Type = -1  ----   perron [G(jw)]
%     Type = 0   ----   freq[G(jw)] 
%     Type = 1   ----   G(jw) 
%     Type = 2   ----   inv(G(jw))
%     Type = 3   ----   I + G(jw)
%     Type = 4   ----   I + inv(G(jw)) 
%
% Vector W contains the frequencies at which the frequency response
% is to be evaluated. The SV matrix has rows which correspond to the 
% singular values in decending order.
%

% R. Y. Chiang & M. G. Safonov 5/16/85
% Copyright (c) 1988 by the MathWorks, Inc.
% All Rights Reserved.
% -------------------------------------------------------------------
%
disp('  ')
disp('          ..... Working ...... please wait .....')
if Type~=1 & Type ~=2 & Type ~=3

[mg] = freqrc(a,b,c,d,w);
[rmg,cmg] = size(mg);
[rb,cb] = size(b);
[rc,cc] = size(c);
gg = ones(rc,cb);
  if (Type == 0)
    sv = mg;
  else
    for is = 1 : cmg
      gg(:) = mg(:,is);
        if (Type == -1)
          sv(:,is) = perron(gg);
        end

        if (Type == 1)
          sv(:,is) = svd(gg);
        end

        if (Type == 2)
          sv(:,is) = svd(inv(gg));
        end

        if (Type == 3)
          sv(:,is) = svd(eye(cb) + gg);
        end

        if (Type == 4)
          sv(:,is) = svd(eye(cb) + inv(gg));
        end
     end
end

else
   if Type ==1
     [sv] = sigma(a,b,c,d,w);
   elseif Type ==2
     [sv] = sigma(a,b,c,d,w,'inv');
   elseif Type==3
     [sv] = sigma(a,b,c,d+eye(size(d)),w);
   else
     [sv]=[];
   end
end


%
% ----- End of SIGMA.M ---- RYC/MGS 5/16/85 %
