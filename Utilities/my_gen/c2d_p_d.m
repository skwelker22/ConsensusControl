function [ad,bd,cd,dd]=c2d_p_d(a,b,c,d,delay,Ts,con_s,trans);
%function [ad,bd,cd,dd]=c2d_p_d(a,b,c,d,delay,Ts,con_s,trans);
% Conversion of a continuous time plant to discrete time with sampling
%   time Ts and delay. the discretization uses Tustin transformation.
%   The delay is added in integer multiples of Ts.
% If con_s = 1, the result is converted back to continuous time
% (default = 0, no conversion)
% trans is a string defining the type of bilinear transformation
%   Tustin (default), Fwdrec, Bwdrec 

if nargin<7;con_s=0;end
if nargin<8;trans='Tustin';end

nT= round(delay/Ts);
if isempty(a)
   ad=[];bd=[];cd=[];dd=d;
else 
   [ad,bd,cd,dd]= bilin(a,b,c,d,1,trans,Ts);
end
[no,ni]=size(d);

for i=1:nT, 
   [ad,bd,cd,dd]=series(zeros(ni,ni),eye(ni,ni),eye(ni,ni),zeros(ni,ni),...
      ad,bd,cd,dd); 
end
if con_s==1
   [ad,bd,cd,dd]= bilin(ad,bd,cd,dd,-1,trans,Ts);
end

