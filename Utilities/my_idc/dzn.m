function dzn=dzn(x,dzh,dzl);
% usage dzn = dzn(x,dzh,[dzl]);
% computes the output of a continuous dead-zone function 
% with input x and thresholds dzh [dzl]. (can be vectors)
% default: dzl=-dzh

if nargin < 3 
dzl=-dzh;
end

dzn=x-sat(x,dzh,dzl);
