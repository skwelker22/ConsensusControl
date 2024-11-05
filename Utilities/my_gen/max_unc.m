function g = optgap(P);

G=minreal(ncf(P));
Wc=gram(G,'c');
Wo=gram(G,'o');
sv=max(eig(Wc*Wo));
sv=min(sv,1);
g=sqrt(1-sv);
