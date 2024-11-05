function [r,rx]=prbs_1(n,a,bl,bt,bias,seed);
% function r=prbs_1(n,a,bl,bt,bias,seed);
% n: PRBS sequence of 2^n points;
% a: minimum time between switches (~1/BW)
% bl: number of leading zeros
% bt: minimum number of trailing zeros
% bias: bias term

if ~isempty(seed)
  rand('seed',seed);
end
tns=2^n;
rx=rand(round((tns-bl-bt)/a-.5),1)-.5+bias;
rx=rx*100000;rx=max(rx,-1);rx=min(rx,1);
rt=kron(rx,ones(a,1));
r=[zeros(bl,1);rt;zeros(tns-bl-length(rt),1)];

