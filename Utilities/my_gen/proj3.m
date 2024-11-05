function proj3=proj3(z,theta,Ra1,Ra2,Ra3,cent1,cent2,cent3,epsst);
%
% USAGE: proj3(z(n),theta(n),Ra-1-2-3-(n x n),cent-1-2-3-(n),epsst)
%
% The function proj3 computes the Lipschitz projection of the
% vectorfield z on three ellipsoids (x-centi)'Rai(x-centi)<=1 with
% boundary layer epsst.
%

epsst2=epsst/2;

thperp1=Ra1*(theta-cent1);
thperp2=Ra2*(theta-cent2);
thperp3=Ra3*(theta-cent3);
dista1=sqrt(abs((theta-cent1)'*thperp1))-1;dista1=max(dista1(1),0);
dista2=sqrt(abs((theta-cent2)'*thperp2))-1;dista2=max(dista2(1),0);
dista3=sqrt(abs((theta-cent3)'*thperp3))-1;dista3=max(dista3(1),0);
dista=sqrt(dista1*dista1+dista2*dista2+dista3*dista3);

thperp=0*(theta);
if dista1>0
thperp=thperp+dista1*thperp1/(dista1+1);
end
if dista2>0
thperp=thperp+dista2*thperp2/(dista2+1);
end
if dista3>0
thperp=thperp+dista3*thperp3/(dista3+1);
end

psi=(dista-epsst2)/epsst2;psi=min(1,max(psi(1),0));
l_tst=psi*z'*thperp;
if l_tst(1) <=0
	proj3=z;
else
	proj3=z-psi*thperp*thperp'*z/(thperp'*thperp);
end
