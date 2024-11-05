function [stchk,MUNCe]=unc_preg(N_P,D_P,cl_spec,frun,AUNC,flag);
%function [stchk,MUNCe]=unc_preg(N_P,D_P,cl_spec,frun,AUNC,flag);
% Produces an estimate of the stability condition for the following
% feedback structure, given estimates of the closed loop sensitivities.   
%                                    |
%                                    V
%       u2---> + --> C_i --> N_P --> + --> (D_P) ----->y1
%              |                                    |
%             -|                                    |
%              |<-----------------------------------|
%
%   N_P,D_P are packed state space descriptions of the plant
%     left factorization
%   cl_spec contains the closed loop specs [T_bw,T_roll,S_bw,S_roll]
%   frun,AUNC(:,2,3) is the coprime factor uncertainty


%  [Adi,Bdi,Cdi,Ddi]=ssinv(A_D,B_D,-C_D,eye(size(D_D))-D_D);

[Adi,Bdi,Cdi,Ddi]=branch(D_P);
[A_N,B_N,C_N,D_N]=branch(N_P);
[ninp,noutp]=size(D_N');nix=min(ninp,noutp);
[A_p,B_p,C_p,D_p]=series(A_N,B_N,C_N,D_N,Adi,Bdi,Cdi,Ddi);

       ssT=awgt_sel(frun,cl_spec(1),cl_spec(2),0);
       ssS=awgt_sel(frun,cl_spec(3),cl_spec(4),-1);

ssP=sv3_5(A_p,B_p,C_p,D_p,1,frun);ssPi=1.0./ssP(nix,:)';
ssN=sv3_5(A_N,B_N,C_N,D_N,1,frun);ssNi=1.0./ssN(nix,:)';
ssDi=sv3_5(Adi,Bdi,Cdi,Ddi,1,frun);ssDi=ssDi(1,:)';
stchk1=ssT.*ssNi.*AUNC(:,2)+ssS.*ssDi.*AUNC(:,3);
stchk=ssT.*ssPi.*ssDi.*AUNC(:,2)+ssS.*ssDi.*AUNC(:,3);

if flag>=0
  loglog(frun,stchk,'y',frun,stchk1,'g');grid
  title('Stability Condition Check (y) (<1)');
end

%------------------------------- Inner loop effective output mult.unc
%sk12=ssk1.*AUNC(:,2);sk12=sk12-sk12.^2+ssk1.*ssk2.*AUNC(:,2).*AUNC(:,3);
%MUNCn=sk12+ssk2.*AUNC(:,2).*sscs./(ssT(nix,:)');
MUNCe=stchk./(1-stchk);
%----------------------------------------------------------------
