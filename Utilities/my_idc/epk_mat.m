function v=epk_mat(v0,i,a);
%function v=epk_mat(v0,i,a);
%replaces the i-th matrix in the model v0 by a

v=v0;
[nn,mm]=size(a);
if i==1
   n=v(2);m=v(3);
   if norm([nn,mm]-[n,m])==0
      v(4:3+n*m)=vector(a);
   else
      disp('error in epk_mat matrix dimensions')
   end
elseif i>1 & i<=v(1)
   ind=2;
   for j=1:i-1
      n=v(ind);m=v(ind+1);ind=ind+2+n*m;
   end
   n=v(ind);m=v(ind+1);ind=ind+2;
   if norm([nn,mm]-[n,m])==0
      v(ind:ind-1+n*m)=vector(a);
   else
	disp('error in epk_mat matrix dimensions');
   end
else
   disp('error index in epk_mat');
end

