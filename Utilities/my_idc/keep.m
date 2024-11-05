function x=keep(x,range)
%function x=keep(fname,range)

eval(['load ',x]);
if (exist('ref')~=1);ref=0*Yp;end
 
t=t(range);t=t-t(1);
Ys=Ys(range,:);
Yp=Yp(range,:);
Up=Up(range,:);
ref=ref(range,:);

   eval(['save ',x,' t Up Ys Yp ref'])

