function G=graph(P,s);

if nargin<2, s='l';end

if s == 'r'
    [n,m]=ncf(P,1);
    G=[n;m];
else
    [n,m] = ncf(P);
    G=[n m];
end
