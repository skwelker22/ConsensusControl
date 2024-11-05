function [pproj3,ds0]=pproj3(theta,Ra1,Ra2,Ra3,cent1,cent2,cent3,epsst);
%
% USAGE: pproj3(theta(n),Ra-1-2-3-(n x n),cent-1-2-3-(n),epsst)
%
% The function pproj3 computes the projection of the parameter
% vector theta on three ellipsoids (x-centi)'Rai(x-centi)<=1 with
% tolerance epsst/2.
%

epsst2=epsst/2;
epsst4=epsst/4;

thperp1=Ra1*(theta-cent1);
thperp2=Ra2*(theta-cent2);
thperp3=Ra3*(theta-cent3);
dista1=sqrt(abs((theta-cent1)'*thperp1))-1;if dista1<0;dista1=0;end
dista2=sqrt(abs((theta-cent2)'*thperp2))-1;if dista2<0;dista2=0;end
dista3=sqrt(abs((theta-cent3)'*thperp3))-1;if dista3<0;dista3=0;end
dista=sqrt(dista1*dista1+dista2*dista2+dista3*dista3);
pproj3=theta;
ds0=dista;

if dista < epsst2
pproj3=theta;
else
	while dista >= epsst2
	pproj3=pproj1(pproj3,Ra1,cent1,epsst4);
	pproj3=pproj1(pproj3,Ra2,cent2,epsst4);
	pproj3=pproj1(pproj3,Ra3,cent3,epsst4);
	thperp1=Ra1*(pproj3-cent1);
	thperp2=Ra2*(pproj3-cent2);
	thperp3=Ra3*(pproj3-cent3);
	dista1=sqrt(abs((pproj3-cent1)'*thperp1))-1;if dista1<0;dista1=0;end
	dista2=sqrt(abs((pproj3-cent2)'*thperp2))-1;if dista2<0;dista2=0;end
	dista3=sqrt(abs((pproj3-cent3)'*thperp3))-1;if dista3<0;dista3=0;end
	dista=sqrt(dista1*dista1+dista2*dista2+dista3*dista3);
	end
end
