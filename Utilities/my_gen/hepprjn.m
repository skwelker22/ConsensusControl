function [proj,d2]=hepprjn(theta,A,B,RA,CE,epsst,distaty,flag);
%
% USAGE: hepprjn(theta,A,B,RA,CE,epsst,distaty,flag)
%
% The function plpprjn computes a projection of the parameter
% vector theta on mh half-spaces A(i,:) * theta <= B(i) 
% and me ellipsoids (x-CEi)'Rai(x-CEi)<=1 (continuous version)
% stopping tolerance epsst (max distance from constraints).
% flag<=1 runs silent, distaty =0/1 selects distance/eqn error
if nargin<7;distaty=1;flag=0;end
if nargin<8;flag=0;end

n=length(theta);mh=length(B);[nx,me]=size(CE);mmh=length(A);mme=length(RA);
if mh==1 & mmh==1
  if A(1)==0 & B(1)==0, mh=0;end
end
if me==1 & mme==1
  if RA(1)==0 & CE(1)==0, me=0;end
end

proj=theta;epsst2=epsst/4;
[d2,dista,D]=hseldist(theta,A,B,RA,CE,distaty,flag);
V=d2(1);
   while dista >= epsst(1)
     for index1=1:mh+me
       if index1<=mh
         if D(index1)>= epsst(1);
           aa=A(index1,:);bb=B(index1);
           proj=proj-aa'*(aa*proj-bb)/(aa*aa');
         end
       else
         if D(index1)>=epsst(1);
           index2=index1-mh;
          proj=pproj1(proj,RA(:,(index2-1)*n+1:index2*n),CE(:,index2),epsst2);
         end
       end
       [d2,dista,D]=hseldist(proj,A,B,RA,CE,distaty,flag);
     end
     Vn=d2(1);
     if dista >epsst(1) & abs(V-Vn) <=epsst(1)*1.e-6
        if flag(1)>0
        disp('Possibly infeasible constraints: terminating projection iterations');
        end
        dista=0;
     end
     V=Vn;  
   end

