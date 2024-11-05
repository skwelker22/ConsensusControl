function [proj,d2,dinf]=hepprjm(theta,A,B,RA,CE,epsst,distaty,flag,f_cont);
%
% USAGE: [proj,d2,dinf]=hepprjm(theta,A,B,RA,CE,epsst,distaty,flag,f_cont);
%
% Computes a projection of the parameter vector theta 
% on mh half-spaces A(i,:) * theta <= B(i) 
% and me ellipsoids (x-CEi)'Rai(x-CEi)<=1. 
%   epsst = stopping tolerance (max distance from constraints; def=1e-5).
%   distaty = 1 for normalized matrix A (norm of each row = 1) 
%   flag <= 1 runs silent, 
%   f_cont = 0 for continuous (in theta) projections (slower)
%   proj = projected point,
%   d2,dinf = error norms; if dinf>epsst 
%             no convergence was achieved (infeasible constraints)

thresh=1.e-12;
if nargin<6;distaty=0;flag=0;epsst=0;f_cont=0;end
if nargin<7;distaty=0;flag=0;f_cont=0;end
if nargin<8;flag=0;f_cont=0;end
if nargin<9;f_cont=0;end
if epsst(1) == 0;epsst=1.e-5;end

n=length(theta);mh=length(B);[nx,me]=size(CE);mmh=length(A);mme=length(RA);
if mh==1 & mmh==1
  if A(1)==0 & B(1)==0, mh=0;end
end
if me==1 & mme==1
  if RA(1)==0 & CE(1)==0, me=0;end
end

if distaty(1) ~=1 & mh >0
   for i=1:mh
     anor=sqrt(A(i,:)*A(i,:)');
       if anor(1) > thresh
         A(i,:)=A(i,:)/anor;B(i)=B(i)/anor;
       else
         A(i,:)=0*A(i,:);B(i)=0;
          if flag(1) >1;disp('ignoring possibly singular constraint');end
       end
   end
end


proj=theta;epsst2=epsst/4;
[d2,dista,D]=hseldist(theta,A,B,RA,CE,1,flag);
V=d2(1);dinf=dista;

   while dista >= epsst(1)

     if f_cont(1)~=0 
       index1=max(find(D==max(D)));
         if index1<=mh
           aa=A(index1,:);bb=B(index1);
           proj=proj-aa'*(aa*proj-bb);
         else
           index2=index1-mh;
          proj=pproj1(proj,RA(:,(index2-1)*n+1:index2*n),CE(:,index2),epsst2);
         end
       [d2,dista,D]=hseldist(proj,A,B,RA,CE,1,flag);
     else
       for index1=1:mh+me
         if index1<=mh
           if D(index1)>= epsst(1);
             aa=A(index1,:);bb=B(index1);
             proj=proj-aa'*(aa*proj-bb);
             [d2,dista,D]=hseldist(proj,A,B,RA,CE,1,flag);
           end
         else
           if D(index1)>=epsst(1);
             index2=index1-mh;
          proj=pproj1(proj,RA(:,(index2-1)*n+1:index2*n),CE(:,index2),epsst2);
             [d2,dista,D]=hseldist(proj,A,B,RA,CE,1,flag);
           end
         end
       end
     end

     Vn=d2(1);dinf=dista;
       if  dista>epsst(1) & abs(V-Vn) <=epsst(1)*1.e-6
          if flag(1) >0
            disp('Possibly infeasible constraints: terminating projection iterations');
          end
          dista=0;
       end
     V=Vn;  
   end

