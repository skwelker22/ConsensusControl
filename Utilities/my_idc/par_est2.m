function [thx,err,WWX]=par_est2(ninp,noutp,filord,FFF,CCC,uf,yf,y,t,cutoff,nc,red_meth);
%function thx=par_est2(ninp,noutp,filord,FFF,CCC,uf,yf,y,t,cutoff,nc,red_meth);
% Called by MIMOID2, for the parameter estimation in system id.

nn=length(uf);

%m2=lsim(-2*del0,1,1,0,1.*u.*u+q2y*y.*y,t);
%m3=ones(m2)+norflag*m2; m2=1+m2;

  if nc==1 | nc==0
    thx=zeros((ninp+2)*max(filord)+ninp*(1-nc),noutp);
  else
    disp('Incorrect spec. of constraints; reset to biproper')
    thx=zeros((ninp+2)*max(filord)+ninp,noutp);
    nc=0;
  end
err=zeros(nn,noutp);

  for i=1:noutp
     if i==1,nftem=0;else,nftem=sum(filord(1:i-1));end
     nftemi=nftem+filord(i);n=filord(i);
     f=FFF(nftem+1:nftemi,nftem+1:nftemi)';
     q=CCC(i,nftem+1:nftemi)';
     cy=eye(n,n); zq=q*0;wu=zeros(nn,ninp*n);
       for ii=1:ninp
         wu(:,(ii-1)*n+1:ii*n)=lsim(f,q,cy,zq,uf(:,ii),t);
       end
     wy=lsim(f,q,cy,zq,yf(:,i),t);
     w0=lsim(f,0*q,cy,zq,0*t,t,q);
      if nc==1
        www=[wu,wy,w0];
      else
        www=[wu,wy,uf,w0];
      end
%   thtem=inv(www'*www)*www'*(y(:,i));
    wwww=www'*www;
    [UW,SW,VW]=svd(wwww);
    tol = length(wwww)*SW(1)*eps;	
    sind = find(diag(SW)>tol);       
    V1=VW(:,sind);S1=diag(1./diag(SW(sind,sind)));
    YY=UW'*www'*(yf(:,i));
    thtem=V1*S1*YY(sind);
    erri=yf(:,i)-www*thtem;norerri=erri'*erri;
      if cutoff(1) > 0
        if red_meth==1
          YY=sqrt(S1)*YY(sind);
          sumerr=0;sind1=length(sind);
            while sumerr <= cutoff(1)*norerri
              sumerr=sumerr+YY(sind1)*YY(sind1);sind1=sind1-1;
            end
          sind1=sind1+1;disp([sind1,length(thtem)]);
          V1=VW(:,1:sind1);S1=S1(1:sind1,1:sind1);
          thtem=V1*sqrt(S1)*YY(1:sind1);
        else
          LW=(0*thtem+cutoff(3));ik=min(i,ninp);
          LW(length(LW)-n+1:length(LW))=...
                   LW(length(LW)-n+1:length(LW))*0+1;
          LW((ik-1)*n+1:ik*n)=LW((ik-1)*n+1:ik*n)*0+1;
          LW(ninp*n+1:ninp*n+n)=LW(ninp*n+1:ninp*n+n)*0+cutoff(2);
            if nc~=1,LW(ninp*n+n+ik)=1;end
          RW=diag(LW);
          Ra=inv(RW)*www'*www*inv(RW)/(max(cutoff(1),1.e-4));
          thtem=orppr1(thtem*0,Ra,RW*thtem,1.e-4*cutoff(1));
          thtem=inv(RW)*thtem;
        end
      end
    thx(:,i)=[thtem;zeros((ninp+2)*(max(filord)-n),1)];
    err(:,i)=y(:,i)-www*thtem;
  end

clg;plot(t,err);pause
    WWX=wwww;
