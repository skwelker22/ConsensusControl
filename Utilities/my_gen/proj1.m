function proj1=proj1(z,theta,Ra,cent,epsst);
%
% USAGE: proj1(z(n),theta(n),Ra(n x n),cent(n),epsst)
%
% The function proj1 computes the Lipschitz projection of the
% vectorfield z on the ellipsoid (x-cent)'Ra(x-cent)<=1 with
% boundary layer epsst.
%

thperp=Ra*(theta-cent);
dista=sqrt(abs((theta-cent)'*thperp))-1;dista=max(dista(1),0);
psi=min(dista/epsst(1),1);
l_tst=psi*z'*thperp;

if l_tst(1) <=0
	proj1=z;
else
	proj1=z-psi*thperp*thperp'*z/(thperp'*thperp);
end
