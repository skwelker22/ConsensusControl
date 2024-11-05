function UNC=unc_bn(fer,fu,fy,ssva,utype,opt_flag);
%function UNC=unc_bn(fer,fu,fy,ssva,utype,opt_flag);
% estimates the uncertainty in sys.id. (used by MIMOID/UNCERF)
% err= estimation error fft, u=input fft (to the uncertainty)
% ssva= singular value weight of an appropriate system
% utype = 1: err/u * ssva, 2: UNC=[UNC1,UNC2], the solution
%       minimizing S1 * UNC1 + S2 * UNC2
%                s.t.  err = UNC1 * u + UNC2 * y
%                where S1 and S2 are the rows of ssva at freq. frun
%         usually, frun is obtained from a preliminary run
%   if no y is supplied, set y=u; if no ssva, set ssva=1;

if nargin<6;opt_flag=0;end
toler=1.e-8;

if utype == 1
   fxx=fer./fu;
   UNC=(fxx).*ssva';
elseif utype == 2
   feu=fer./fu;
   fey=fer./fy;
   fyu=fy./fu;
   ab=((ssva(1,:)./ssva(2,:))').*fyu;
   den_temp=sqrt(1+ab.*ab);
   ind_y=find(ab>1);ind_u=find(ab<=1);
      if opt_flag==0
        UNCu=0*ab+toler;UNCy=UNCu;
        UNCu(ind_u)=feu(ind_u);
        UNCy(ind_y)=fey(ind_y);
      elseif opt_flag==2
        UNCu=feu./den_temp;
        UNCy=fey./den_temp.*ab;
     else
        UNCu=feu./(1+ab);
        UNCy=fey.*ab./(1+ab);
      end
   UNC=[UNCu,UNCy];
end
