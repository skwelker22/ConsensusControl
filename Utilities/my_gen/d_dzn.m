function dd = d_dzn(d,dd,lam,dzn)
% for i=2:length(d);dd(i)=d_dzn(d(i),dd(i-1),10,.1);end;plot(t,d,t,dd)

e=d-dd;
if e>dzn; 
   e=e-dzn;
elseif e<-dzn
   e=e+dzn;
else
   e=0;
end

dd=dd+lam*e^3/(1+lam*e^2);

