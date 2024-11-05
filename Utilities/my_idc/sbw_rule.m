function sbw=sbw_rule(cl_sp,cns_fct);
% function sbw=sbw_rule(cl_sp,cns_fct);
% provides an ad-hoc rule for the selection of the closed-loop
% sensitivity cut-off frequency, based on the desired closed-loop
% bandwidth and roll-off rates that are specified in cl_sp
% INPUTS:  cl_sp = [T_bw, T_rolloff, xxx, S_rolloff, xxx, AGRS]
% OUTPUT:  S_bw

if nargin<2;cns_fct=1;end
if length(cl_sp)>=6;AGRS=cl_sp(6);else;AGRS=1;end
cns_fct=cns_fct*AGRS;

if length(cl_sp)==1
   sbw=cl_sp*.8*cns_fct;
elseif length(cl_sp)>=4
   z1=max(cl_sp(2)-1,0);z2=max(cl_sp(4)-1,0);
   if cl_sp(4)>2;z3=cl_sp(4)*2.2/3;else;z3=cl_sp(4);end
   sbw=cl_sp(1)*(0.9^z1)*(0.65^z2)^(z3)*cns_fct;
   if cl_sp(1)>3;sbw=sbw*.8;end
   if cl_sp(4)>=2;sbw=sbw*.9;end
else
   z1=max(cl_sp(2)-1,0);
   sbw=cl_sp(1)*(0.9^z1)*cns_fct;
end
