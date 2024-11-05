function y=lsim_nan(a,b,c,d,u,t,x0,ic);

% function  y=lsim_nan(a,b,c,d,u,t,x0,ic);
% lsim with discontinuous data

if isobject(a)
    H=ss(a); 
    if nargin >=3, t=c; end
    if nargin >=4, x0=d; else x0=[]; end
    if nargin >=5, ic=u;, else ic=0; end
    u=b;
else
    H=ss(a,b,c,d);
    if nargin < 7;x0=[];end
    if nargin < 8; ic=0;end
end

DC=dcgain(H);
n=length(u);  [noutp,ninp]=size(H.d);   y=zeros(n,noutp);
N = length(H.a);
nan_par=[find(isnan(u(:,1)));n+1];    n_nan=length(nan_par);  
for i=1:n_nan-1;   y(nan_par(i),:)=y(nan_par(i),:)*NaN; end
nan_ind1=1;

if ic ~=0;  if H.Ts == 0; L = [H.a; H.c]; else L=[H.a-eye(size(H.a)); H.c];end;end

for i_nan=1:n_nan
   nan_ind=[nan_ind1:nan_par(i_nan)-1];
   z0=[zeros(N,1)-H.b*(u(nan_ind1,:)'); DC*u(nan_ind1,:)'-H.d*(u(nan_ind1,:)')];
   nan_ind1=nan_par(i_nan)+1;
   X0=[];
   if ic ~= 0; X0 = L\z0; X0=X0(1:N); end
   if ~isempty(x0), X0=x0(:,i_nan);end
   if isempty(X0) 
       if H.Ts == 0,
        y(nan_ind,:)=lsim(H,u(nan_ind,:),t(nan_ind)-t(nan_ind(1)));
       else
        y(nan_ind,:)=lsim(H,u(nan_ind,:));
       end           
   else
       if H.Ts == 0
        y(nan_ind,:)=lsim(H,u(nan_ind,:),t(nan_ind)-t(nan_ind(1)),X0);
       else
           y(nan_ind,:)=lsim(H,u(nan_ind,:),[],X0);
       end
   end
end
