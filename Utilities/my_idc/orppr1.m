function orppr1=orppr1(theta,Ra,cent,epsst,Niter);
%
% USAGE: orppr1(theta(n),Ra(n x n),cent(n),epsst,Niter)
%
% The function orppr1 computes the orthogonal projection of the parameter
% vector theta on the ellipsoid (x-cent)'Ra(x-cent)<=1 with
% tolerance epsst.
%
if nargin<5;Niter=100;end
kiter=0;

NN=length(theta);INN=eye(NN,NN);

epsst2=1+epsst;
xi=theta-cent;
AA=(xi'*Ra*xi);
if AA < epsst2
orppr1=theta;
else
b=sqrt(AA);
Z=(INN+b*Ra);
z1=Z\xi;r1=Ra*z1;z2=Z\r1;
f=xi'*z2-1;
	while abs(f) >= epsst
	kiter=kiter+1;
      r2=Ra*z2;
	grad=-2*xi'*(Z\r2);
	db=-f/grad;
		while (b+db) <= 0
		db=db/2;
		end
	b=b+db;
	Z=(INN+b*Ra);
      z1=Z\xi;r1=Ra*z1;z2=Z\r1;
	f=xi'*z2-1;
        if kiter>Niter;f=0;end
	end
orppr1=z1+cent;
end
%kiter
