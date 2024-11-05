function a=unpk_mat(v,i_mat);
% function a=unpk_mat(v,i_mat);

n_mat = v(1);
if i_mat > n_mat;
   disp('error index in unpk_mat');
   a=[]; return
end

ind=2;
for i=1:i_mat
   n=v(ind);m=v(ind+1);ind=ind+2;
   a=unvector(v(ind:ind+n*m-1),n,m);
   ind=ind+n*m;
end
  