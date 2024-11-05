function [J,DJ,DDJ]=pe_obj(x,indc,cutoff,f,q,cy,zq,RWi,N_A,t,u,y,ww);
%  function [J,DJ,DDJ]=proj_obj(x,indc);

J=0;DJ=0;DDJ=0;

[n,nx]=size(f);[ndat,ninp]=size(u);
R=RWi*N_A;sigma=1e-2*cutoff(1);
th=R*x;
thd=th(1:ninp);
thN=unvector(th(1+ninp:ninp*n+ninp),n,ninp);
thD=th(ninp*n+ninp+1:ninp*n+ninp+n);
th0=th(ninp*n+ninp+n+1:ninp*n+ninp+n+n);

yh=lsim_nan(f'+thD*q',thN+thD*thd',q',thd',u,t,th0);
eh=yh-y;

J=(eh'*eh)/2+sigma*th'*th/2;
if indc(1)==1
   ww(:,ninp*(1+n)+1:ninp*(1+n)+n)=lsim_nan(f,q,cy,zq,yh,t(:)-t(1));
   DJ=(eh'*ww*R)'+sigma*(th'*R)';
   DDJ=R'*ww'*ww*R+sigma*R'*R;
end
