function x=cmp_mdls(m1,m2,m3,m4,m5);

f1=[];f2=[];f3=[];f4=[];f5=[];
s1=[];s2=[];s3=[];s4=[];s5=[];
eval([' load ',m1]);
f1=frw;s1=sspla;
if nargin>=2
  eval([' load ',m2]);f2=frw;s2=sspla;
end
if nargin>=3
  eval([' load ',m3]);f3=frw;s3=sspla;
end
if nargin>=4
  eval([' load ',m4]);f4=frw;s4=sspla;
end
if nargin>=5
  eval([' load ',m5]);f5=frw;s5=sspla;
end
loglog(f1,s1,'y',f2,s2,'r',f3,s3,'b',f4,s4,'g',f5,s5,'c')
title('SV sequence: y,r,b,g,c')
