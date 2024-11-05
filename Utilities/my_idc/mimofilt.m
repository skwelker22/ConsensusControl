function y=mimofilt(num,den,u,t,ic_pts);
% function y=mimofilt(num,den,u,t,ic_pts);
% multivariable filtering of each component of u
%     with the transfer function num/den.
%     (t,u) is the input function where t is
%       equally spaced time samples and u contains
%       the corresponding input values in columns
%     ic_pts: optional; sets initial conditions 
%       to the steady-state of the average of the first
%       ic_pts points (def=1)

temp=version;
if str2num(temp(1))>=6;
    if isobject(num);[num,den]=tfdata(num,'v');end
end

if nargin < 5; ic_pts=1;end
if isempty(ic_pts); ic_pts=1;end
x0=zeros(length(den)-1,1);
[n,m]=size(u);
DC=num(length(num))/den(length(den));
y=u*DC;

[a,b,c,d]=tf2ss(num,den);
[rb,cb]=size(b);ry=1;
A=[-a;c];
B=[b,zeros(rb,ry);-d,eye(ry,ry)];
L=pinv(A)*B;
CCA=[c;c*a];

nan_par=[find(isnan(u(:,1)));n+1];
n_nan=length(nan_par);nan_ind1=1;
for i_nan=1:n_nan
   nan_ind=[nan_ind1:nan_par(i_nan)-1];
   for i=1:m
      if ic_pts >0
         ym=mean(y(nan_ind1:nan_ind1+ic_pts,i));
         um=mean(u(nan_ind1:nan_ind1+ic_pts,i));
         x0=L*[um;ym];
         %    x0=c\(ym-d*um);
      elseif ic_pts<-1
         t1=t(nan_ind1:nan_ind1-ic_pts);t1=t1-t1(1);
         ab=pinv([0*t1+1,t1])*u(nan_ind1:nan_ind1-ic_pts,i);
         bb=[(DC-d)*ab(1)-c*inv(a*a)*b*ab(2);(DC-d)*ab(2)-c*b*ab(1)];
         x0=pinv(CCA)*bb;
      end
      y(nan_ind,i)=lsim4(a,b,c,d,u(nan_ind,i),t(nan_ind)-t(nan_ind(1)),x0);
   end
   nan_ind1=nan_par(i_nan)+1;
end
