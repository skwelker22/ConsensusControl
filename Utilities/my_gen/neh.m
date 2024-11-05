function [as,bs,cs,ds,err]= neh(a,b,c,d)
% function [as,bs,cs,ds,err]= neh(a,b,c,d)
% Nehari approximation 

[ra,ca]= size(a);
[ee,dd]= eig(a);
indstab = find(real(dd) < 0);
m = length(indstab);

if m == ra  % ------ If completely stable :
  as=a; bs=b; cs=c; ds=d; err=0;
elseif m==0 % ------ If completely unstable :
  [al,bl,cl,dl,ar,br,cr,dr,aug] = ohkapp(-a,-b,c,d,1,0);
  as= -ar; bs= -br; cs= cr; ds= dr; err= aug(1,1);
elseif m > 0 & m < ra  % ------ If having both stable & unstable parts :
  [al,bl,cl,dl,ar,br,cr,dr] = stabproj(a,b,c,d); 
  [arl,brl,crl,drl,ar,br,cr,dr,aug] = ohkapp(-ar,-br,cr,dr,1,0);
  [as,bs,cs,ds]= addss(al,bl,cl,dl,-ar,-br,cr,dr);
  err= aug(1,1);
end
% svout = sigma(ss(addss(a,b,c,d,as,-bs,cs,-ds)));
