function [thx,err,nan_ind]=par_esta(ninp,noutp,filord,FFF,CCC,uf,yf,t,cutoff,nc,red_meth);
%function thx=par_est(ninp,noutp,filord,FFF,CCC,uf,yf,t,cutoff,nc,red_meth);
% Called by MIMOID, for the parameter estimation in system id.

tolopt=1e-4;Nmax=50;
COV=[];COV_R=[];COV_cns=[];

nn=length(uf);
nan_ind=[find(isnan(uf(:,1)));nn+1];
n_nan=length(nan_ind)-1;
n_tot=nn-n_nan;
i_nnan=find(~isnan(uf(:,1)));

ncx=min(1,nc);
thx=zeros((ninp+2+n_nan)*max(filord)+ninp,noutp);
err=zeros(nn,noutp)*NaN;

%m2=lsim4(-2*del0,1,1,0,1.*u.*u+q2y*y.*y,t);
%m3=ones(m2)+norflag*m2; m2=1+m2;

  for i=1:noutp
disp(['Estimating models for output ',num2str(i)]) 
     if i==1,nftem=0;else,nftem=sum(filord(1:i-1));end
     nftemi=nftem+filord(i);n=filord(i);
     f=FFF(nftem+1:nftemi,nftem+1:nftemi)';
     q=CCC(i,nftem+1:nftemi)';
     cy=eye(n,n); zq=q*0;
     n_th=(ninp+2+n_nan)*n+ninp;
     www=zeros(n_tot,n_th);
     YF=zeros(n_tot,1);
     i_nan1=1;

    thtem=zeros(n_th,1);
    LW=(thtem+cutoff(3));   ik=min(i,ninp);
    LW(ik)=1;
    i_xn=ninp+(ik-1)*n+1:ninp+ik*n;      LW(i_xn)=LW(i_xn)*0+1*cutoff(5);
    i_xp=ninp*(1+n)+1:ninp*(1+n)+n;      LW(i_xp)=LW(i_xp)*0+cutoff(2);
    i_x0=n_th-n*(n_nan+1)+1:n_th;        LW(i_x0)=LW(i_x0)*0+1;
    RW=diag(LW);RWi=diag(1./LW);

     if nc>=1 
        Acns=zeros((nc)*ninp,n_th);
	  Acns(1:ninp,1:ninp)=eye(ninp,ninp);
        if nc>1
        for jnc = 1:nc-1
           qq=((f^(jnc-1))*q)';
           for jcns = 1:ninp
              irow=ninp+jcns+(jnc-1)*ninp;
              icol=ninp+1+n*(jcns-1):ninp+n*jcns;
              Acns(irow,icol)=qq;
           end
        end
        end
        N_Acns=null(Acns*RWi);
     else 
        N_Acns=1;
     end

     for i_nan=1:n_nan+1
       ind_nan=[i_nan1:nan_ind(i_nan)-1];
       ind_nan2=ind_nan-i_nan+1;
       i_nan1=nan_ind(i_nan)+1;
       YF(ind_nan2)=yf(ind_nan,i);

       www(ind_nan2,1:ninp)=uf(ind_nan,:);
       for ii=1:ninp
         www(ind_nan2,(ii-1)*n+1+ninp:ninp+ii*n)=lsim_nan(f,q,cy,zq,uf(ind_nan,ii),t(ind_nan)-t(ind_nan(1)));
       end
       www(ind_nan2,ninp*(1+n)+1:ninp*(1+n)+n)=lsim_nan(f,q,cy,zq,yf(ind_nan,i),t(ind_nan)-t(ind_nan(1)));
       w0=lsim_nan(f,0*q,cy,zq,0*t(ind_nan),t(ind_nan)-t(ind_nan(1)),q);
       nan_loc=ninp+(i_nan-1)*n;
       www(ind_nan2,1+(ninp+1)*n+nan_loc:(ninp+1)*n+nan_loc+n)=w0;
     end

%   thtem=inv(www'*www)*www'*(yf(:,i));

    wwww=www'*www;
    wwww=N_Acns'*RWi'*wwww*RWi*N_Acns;
    [UW,SW,VW]=svd(wwww);
    tol = length(wwww)*SW(1)*eps;	
    sind = find(diag(SW)>tol);       
    V1=VW(:,sind);S1=diag(1./diag(SW(sind,sind)));
    YY=UW'*N_Acns'*RWi'*(www'*YF);
    thtem=V1*S1*YY(sind);
    ZLS=www*RWi*N_Acns*thtem;Z1=wwww*thtem;
    erri=YF-ZLS;norerri=erri'*erri;

      if cutoff(1) > 0
        if red_meth==1
          YY=sqrt(S1)*YY(sind);
          sumerr=0;sind1=length(sind);
            while sumerr <= cutoff(1)*norerri;
              sumerr=sumerr+YY(sind1)*YY(sind1);sind1=sind1-1;
            end
          sind1=sind1+1;disp([sind1,length(thtem)])
          V1=V1(:,1:sind1);
          thtem=V1*(sqrt(S1(1:sind1,1:sind1))*YY(1:sind1));
        else
          Ra=wwww/(max(cutoff(1)*norerri,1.e-4));
          thtem=orppr1(thtem*0,Ra,thtem,1.e-4*cutoff(1)*norerri);
        end
      end
    thtem=pe_newt(thtem,tolopt,Nmax,cutoff,f,q,cy,zq,RWi,N_Acns,t,uf,yf(:,i),www);
    thtem=RWi*N_Acns*thtem;

    thx(:,i)=[thtem;zeros((ninp+2+n_nan)*(max(filord)-n),1)];
    err(i_nnan,i)=YF-www*thtem;
    COV=pk_mat(COV,wwww);COV_R=pk_mat(COV_R,RWi);COV_cns=pk_mat(COV_cns,N_Acns);
  end
