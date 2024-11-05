function y=mimofilt6(sys,u,t,ic_pts);
% function y=mimofilt6(sys,u,t,ic_pts);
% multivariable filtering of each component of u
%     with the transfer function num/den.
%     (t,u) is the input function where t is
%       equally spaced time samples and u contains
%       the corresponding input values in columns
%     ic_pts: optional; sets initial conditions 
%       to the steady-state of the average of the first
%       ic_pts points (def=1)


sys=ss(sys); [nout,nin]=size(sys.d);
x0=zeros(length(sys.a),1);
[n,m]=size(u);
DC=sys.d-sys.c*pinv(sys.a)*sys.b;
y=u*DC';

[rb,cb]=size(sys.b);ry=1;
A=[-sys.a;sys.c];
B=[sys.b,zeros(rb,nout);-sys.d,eye(nout,nout)];
L=pinv(A)*B;
CCA=[sys.c;sys.c*sys.a];
nan_par=[find(isnan(u(:,1)));n+1];
n_nan=length(nan_par);nan_ind1=1;
for i_nan=1:n_nan
   nan_ind=[nan_ind1:nan_par(i_nan)-1];
   ym=mean(y(nan_ind1:nan_ind1+ic_pts,:));
   um=mean(u(nan_ind1:nan_ind1+ic_pts,:));
   x0=L*[um';ym'];
   y(nan_ind,:)=lsim4(sys.a,sys.b,sys.c,sys.d,u(nan_ind,:),t(nan_ind)-t(nan_ind(1)),x0);

   nan_ind1=nan_par(i_nan)+1;
end
