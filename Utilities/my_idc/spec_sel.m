function cl_spec=spec_sel(frun,MUNC,S_roll,R_A,N_P,D_P,AUNC,flag,TR);
% fcn cl_spec=spec_sel(frun,MUNC,S_roll,R_A,N_P,D_P,AUNC,flag,TR);
% to determine closed-loop specs with coprime factor uncertainty
% computations.  Checks for user input if flag>=1

if length(R_A)==1;R_A=[R_A,1];end
cl_spec=bw_app(frun,MUNC,S_roll,R_A(1),R_A(2));

stchk=unc_preg(N_P,D_P,cl_spec,frun,AUNC,flag);
nfh=max(find(frun < 8*cl_spec(1)));
st_con=max(stchk(1:nfh));
Kiter=0;
if st_con < TR(1); while st_con<TR(1)
   Kiter=Kiter+1;
   cl_spec(1)=1.1*cl_spec(1);
   cl_spec(3)=sbw_rule(cl_spec);
   stchk=unc_preg(N_P,D_P,cl_spec,frun,AUNC,-1);
   st_con=max(stchk(1:nfh));
   if Kiter>=4;st_con=TR(1)+1;end
end;end
Kiter=0;st_con=max(stchk(1:nfh));
if st_con > TR(2); while st_con>TR(2)
   Kiter=Kiter+1;
   cl_spec(1)=0.9*cl_spec(1);
   cl_spec(3)=sbw_rule(cl_spec);
   stchk=unc_preg(N_P,D_P,cl_spec,frun,AUNC,-1);
   st_con=max(stchk(1:nfh));
   if Kiter>=5;st_con=TR(2)-1;end
end;end
   stchk=unc_preg(N_P,D_P,cl_spec,frun,AUNC,0);

disp('Spec_sel: closed loop specs [T_bw,T_roll,S_bw,S_roll,R,A]')
disp(cl_spec)
s_done=0;

if flag>=1
while s_done==0
     bwdth=input('Spec_sel: Enter target loop T,S specs:  ');
     if length(bwdth)>0;
        cl_spec(1:length(bwdth))=bwdth;
        if cl_spec(3)==0;cl_spec(3)=sbw_rule(cl_spec);end
        stchk=unc_preg(N_P,D_P,cl_spec,frun,AUNC,1);
        s_done=input('Spec_sel: Done with target loop T,S specs? (0=no) [1]   ');
     else
        s_done=1;
     end

   R_A=inv_spec(cl_spec,frun,MUNC);
     if length(cl_spec)<6
       cl_spec(5:6)=R_A;
     else
       if abs(log(R_A(1)/cl_spec(5)))+abs(R_A(2)-cl_spec(6))>0.002
         cl_spec(5:6)=R_A;
       end
     end

end
end

