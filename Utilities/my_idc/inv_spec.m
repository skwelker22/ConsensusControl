function R_A=inv_spec(cl_spec,frun,MUNC);
% function function R_A=inv_spec(cl_spec,frun,MUNC);
% inverts quick spec to find risk and agrs factors


   BW0=cl_spec(1);
   if length(cl_spec)>=5;R=cl_spec(5);else;R=1;end
   temp_s=bw_app(frun,MUNC,cl_spec(4),R,1,-cl_spec(2),0);
   BW1=temp_s(1);
   err=log(BW1/BW0);
kiter=0;
while abs(err)>.001
   kiter=kiter+1;
   R=exp(log(R)-err/(1+abs(err)));
   temp_s=bw_app(frun,MUNC,cl_spec(4),R,1,-cl_spec(2),0);
   BW1=temp_s(1);
   err=log(BW1/BW0);
   if kiter > 50;err=0;disp('inv_spec: max. iterations exceeded');end
end


temp_s=sbw_rule([cl_spec(1:4),1,1]);
AGRS=cl_spec(3)/temp_s;
R_A=[R,AGRS];
