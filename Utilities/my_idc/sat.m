function sat=sat(x,sh,sl);
% usage: sat = sat(x,sh,[sl]);
% computes the output of a saturation function 
% with input x and thresholds sh, [sl]. (Can be vectors but not matrices)
% (default: sl=-sh)

if nargin < 3
sl=-sh;
end

sat=x;[n,m]=size(x);nmax=max(n,m);
shh=0*x+sh;sll=0*x+sl;
z1=x-shh;ind=find(z1>0);
sat(ind)=shh(ind);
z1=x-sll;ind=find(z1<0);
sat(ind)=sll(ind);

inan=find(isnan(x));sat(inan)=nan;

