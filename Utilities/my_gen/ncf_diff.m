function [D_S,N_S]=ncf_diff(N_1,N_2,D_1,D_2,SW_N,red_tol);

%function [D_S,N_S]=ncf_diff(N_1,N_2,D_1,D_2,SW_N,red_tol);

[An,Bn,Cn,Dn]=branch(N_1);
if ~isempty(SW_N)
   [As,Bs,Cs,Ds]=branch(SW_N);
   [An,Bn,Cn,Dn]=series(As,Bs,Cs,Ds,An,Bn,Cn,Dn);
end
S_n=mksys(An,Bn,Cn,Dn,'ss');
[S_i,S_p,S_o]=iofc(S_n);S_oi=ssinv(S_o);
[Ani,Bni,Cni,Dni]=w_sysred(S_i,[],[],0,0,[3,red_tol]);
[An2,Bn2,Cn2,Dn2]=branch(N_2);
if ~isempty(SW_N)
   [Ax,Bx,Cx,Dx]=series(As,Bs,Cs,Ds,An2,Bn2,Cn2,Dn2);
   [Ax,Bx,Cx,Dx]=series(-Ani',-Cni',Bni',Dni',Ax,Bx,Cx,Dx);
else
   [Ax,Bx,Cx,Dx]=series(-Ani',-Cni',Bni',Dni',An2,Bn2,Cn2,Dn2);
end

disp('ncf-diff, Solving Nehari...')
[AQ,BQ,CQ,DQ,sigr]=nehari(Ax,Bx,Cx,Dx);
disp(['ncf-diff, Nehari sigma = ',num2str(sigr)])
[Ai,Bi,Ci,Di]=branch(S_oi);
[AQ,BQ,CQ,DQ]=series(Ai,Bi,Ci,Di,AQ,BQ,CQ,DQ);

[A1,B1,C1,D1]=branch(D_1);[A2,B2,C2,D2]=branch(D_2);
[AD,BD,CD,DD]=series(A1,B1,C1,D1,AQ,BQ,CQ,DQ);
[AD,BD,CD,DD]=addss(AD,BD,CD,DD,A2,-B2,C2,-D2);
STEM=mksys(AD,BD,CD,DD,'ss');[AD,BD,CD,DD]=w_sysred(STEM,[],[],0,0,[3,red_tol]);
D_S=mksys(AD,BD,CD,DD,'ss');

[A1,B1,C1,D1]=branch(N_1);[A2,B2,C2,D2]=branch(N_2);
[AD,BD,CD,DD]=series(A1,B1,C1,D1,AQ,BQ,CQ,DQ);
[AD,BD,CD,DD]=addss(AD,BD,CD,DD,A2,-B2,C2,-D2);
STEM=mksys(AD,BD,CD,DD,'ss');[AD,BD,CD,DD]=w_sysred(STEM,[],[],0,0,[3,red_tol]);
N_S=mksys(AD,BD,CD,DD,'ss');
