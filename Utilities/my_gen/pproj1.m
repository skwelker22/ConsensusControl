function [pproj1,dist0]=pproj1(theta,Ra,cent,epsst);
%
% USAGE: pproj1(theta(n),Ra(n x n),cent(n),epsst)
%
% The function pproj1 computes the projection of the parameter
% vector theta on the ellipsoid (x-cent)'Ra(x-cent)<=1 with
% tolerance epsst/2.
%

epsst2=1+epsst/2;
xi=theta;
xiperp=Ra*(xi-cent);
AA=sqrt(abs((xi-cent)'*xiperp));
dist0=AA-1;if dist0(1)<0;dist0=0;end
	while AA(1) >= epsst2(1)
	xiperp=Ra*(xi-cent);
	AA=sqrt(abs((xi-cent)'*xiperp));
	xi=xi-(1-1/AA)*xiperp*xiperp'*(xi-cent)/(xiperp'*xiperp);
	end
pproj1=xi;
