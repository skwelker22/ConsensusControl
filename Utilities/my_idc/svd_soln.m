function [thx,err,WWX]=svd_soln(www,yf,cutoff);
%function [thx,err,WWX]=svd_soln(www,yf,cutoff);
% Called by MIMOID, for the parameter estimation in system id.

    wwww=www'*www;
    [UW,SW,VW]=svd(wwww);
    tol = length(wwww)*SW(1)*eps;	
    sind = find(diag(SW)>tol);       
    V1=VW(:,sind);S1=diag(1./diag(SW(sind,sind)));
    YY=UW'*www'*(yf);
    thtem=V1*S1*YY(sind);
    erri=yf-www*thtem;norerri=erri'*erri;
      if cutoff(1) > 0
          YY=sqrt(S1)*YY(sind);
          sumerr=0;sind1=length(sind);
            while sumerr <= cutoff(1)*norerri
              sumerr=sumerr+YY(sind1)*YY(sind1);sind1=sind1-1;
            end
          sind1=sind1+1;disp([sind1,length(thtem)]);
          V1=VW(:,1:sind1);S1=S1(1:sind1,1:sind1);
          thtem=V1*sqrt(S1)*YY(1:sind1);
      end
    thx=thtem;
    err=yf-www*thtem;
    WWX=wwww;
