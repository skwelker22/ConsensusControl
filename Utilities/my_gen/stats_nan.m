function [M,S,N]=stats_nan(U);

[n,m]=size(U);
N=zeros(1,m); M=N*nan; S=M;

for i=1:m
    r=U(:,i);
    ind=find(~isnan(r));
    N(i)=length(ind);
    if (N(i) >0)
        M(i)=mean(r(ind));
        S(i)=std(r(ind));
    end
end
