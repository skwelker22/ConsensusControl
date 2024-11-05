function [ce,Ra,Rainv,inx]=setmem(ce,Ra,Rainv,ym,z,nbd,epsi,inx);
%
% USAGE: 
% [ce(n),Ra(n x n),Rainv(n x n),inx]=
%              setmem(ce(n),Ra(n x n),Rainv(n x n),
%                     ym(1),z(n),nbd(1),epsi(1),inx)
%
% The function setmem computes the updated bounding ellipsoid given
% the initial ellipsoid      (x-ce)'Rainv(x-ce)<=1, 
% the measurement            ym = theta_*'z + noise,
% with an absolute noise bound = nbd.
% epsi : the ellipsoid is updated if  det(Ra_new)<=(1-epsi)det(Ra_old).
% On return, inx=1 if an update occured, inx=-1 if the constraints are
%            inconsistent, and inx=0 otherwise.
% 
if inx == -1
    disp('WARNING (setmem): inconsistent constraints')
else
    inx=0;
    NT=length(ce);
    errx=(ym-ce'*z);
    g=z'*Ra*z;
    beta3=NT*(nbd*nbd-errx*errx)-g;
    if beta3<0
        beta1=(NT-1)*g*g*nbd*nbd;
        beta2=g*((2*NT-1)*nbd*nbd+errx*errx-g);
        discr=beta2*beta2-4*beta1*beta3;
        if discr > 0
            q=(-beta2+sqrt(discr))/2/beta1;if q<0;q=0;end
        else
            q=0;
        end
        aa=1+q*(nbd*nbd-errx*errx)+q*q*errx*errx*g/(1+q*g);
        if aa < 0
            q=0;
            inx=-1;
            disp('WARNING (setmem): inconsistent constraints')
        else
            if aa^NT/(1+q*g) > 1-epsi
                q=0;
            else
                inx=1;
            end
        end
        aa=1+q*(nbd*nbd-errx*errx)+q*q*errx*errx*g/(1+q*g);
        ce=ce+Ra*q*errx*z/(1+q*g);
        Ra=aa*Ra-q*aa*Ra*z*z'*Ra/(1+q*g);
        Rainv=(Rainv+q*z*z')/aa;
    end
end
