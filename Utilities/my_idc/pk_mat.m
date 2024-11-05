function v=pk_mat(v0,a);
%function v=app_mdl(v0,a);
% appends the matrix a to the model v0

[n,m]=size(a);
if isempty(v0)
   v=[0;n;m;vector(a)];
else
   v=[v0;n;m;vector(a)];
end
v(1)=v(1)+1;
