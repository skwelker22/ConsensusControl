% Consensus of a Class of Second-Order Multi-Agent Systems
% With Time-Delay and Jointly-Connected Topologies
% Peng Lin and Yingmin Jia
%% startup and format
clear all; close all; home;

%format
set(0, 'DefaultTextFontName', 'Arial', 'DefaultTextFontSize', 20, 'DefaultTextFontWEight', 'Bold');
set(0, 'DefaultAxesFontName', 'Arial', 'DefaultAxesFontSize', 20, 'DefaultAxesFontWeight', 'Bold', 'DefaultAxesLineWidth', 1.5);
set(0, 'DefaultLineLineWidth', 3, 'DefaultLineMarkerSize', 10);
format long g;

%add yalmip and optimizer settings
addpath(genpath('../YALMIP-master'));
addpath('C:\Program Files\Mosek\10.2\toolbox\r2017a');
ops=sdpsettings('solver','mosek','verbose',0); optTol=1e-6;

%% start simulation
%settings from example
k1=2; nAgents=6; aij=0.5;
A=[-k1/2, k1/2; k1/2, -k1/2];
B=[0, 0; 2/k1, 0];

%lmi variables
tau0=1e-3; tauF=2; tauStep=0.5e-1; nTauSamps=ceil((tauF-tau0)/tauStep);
tauVec=linspace(tau0,tauF,nTauSamps); lenTau=length(tauVec);
sdpvar gamma

%loop over all possible switching conditions
nTop=4; [Qi_top, gamma_i_top, convergenceTop]=deal(cell(nTop,1));

%set tau cell
allTauConvergence=cell(lenTau,1);

for iTau=1:lenTau
    %get current tau
    tau=tauVec(iTau);

    %for each network configuration, check stability
    for ss=1:nTop
        %start at topology a
        [L_thetaBar,L_theta,d_sigI,nComponents]=...
                                    GetGraphTopology(aij, ss);

        %setup LMI for the double-integrator system
        %solve for each connected component
        %init optimization vectors to save off data
        convergenceFlags=nan(nComponents,1);
        Qi=cell(nComponents,1); gamma_i=zeros(nComponents,1);

        %init L_sig_i index
        dd=1;
        for ii=1:nComponents
            %current component
            d_sig=d_sigI(ii); two_dSig=2*d_sig;
            Id_sig=eye(d_sig);

            %for each iteration, form Q_sigI
            Q_sigI=sdpvar(d_sig,d_sig);

            %form phi_sigI
            lSig0=dd; lSigf=dd+d_sig-1;

            %increment dSigIx
            dd=dd+d_sig;
            L_sigI=L_thetaBar(lSig0:lSigf,lSig0:lSigf);

            %if the sub laplacian is all zeros, break
            if(all(L_sigI(:)==0))
                break;
            end

            %form sub term
            phi_sigI=kron(Id_sig,A)-kron(L_sigI,B);

            %form lmi terms
            psi_sigI11=gamma*(phi_sigI+phi_sigI')+tau*(kron(Id_sig,A*A));
            psi_sigI12=kron(Q_sigI,[1;0])-tau*kron(L_sigI,[1;-1]); psi_sigI21=psi_sigI12';
            psi_sigI13=gamma*kron(L_sigI,B); psi_sigI31=psi_sigI13';
            psi_sigI22=(4*tau/k1^2)*L_sigI*L_sigI-2*Q_sigI;
            psi_sigI23=kron(-Q_sigI, [1;0]'); psi_sigI32=psi_sigI23';
            psi_sigI33=-eye(two_dSig)/tau;

            %form lmi
            PSI_SIG=[psi_sigI11, psi_sigI12, psi_sigI13;
                psi_sigI21, psi_sigI22, psi_sigI23;
                psi_sigI31, psi_sigI32, psi_sigI33];

            %form Hi
            [W_sig,~]=eig(L_sigI);
            Hi_barT=blkdiag(kron(W_sig', eye(2)), W_sig', kron(W_sig',eye(2)));
            Hi_bar=blkdiag(kron(W_sig, eye(2)), W_sig, kron(W_sig,eye(2)));

            %form constraints
            %note: the rank constraint isn't gaurunteed using the nuclear norm
            % in order to actually enforce the rank constraint a non-convex
            % solver should be used.
            Constraints=[(Hi_barT * PSI_SIG * Hi_bar)<=optTol, ...
                Q_sigI>=optTol];
            %Q_sigI*ones(d_sig,1)<=optTol.*ones(d_sig,1)];
            %norm(Q_sigI,'nuclear')<=(d_sig-1), ...


            %objective
            Objective=[];

            %solve optimization
            optOut=optimize(Constraints, Objective, ops);

            %keep track of gamma, Qi, and feasible solutions
            convergenceFlags(ii)=optOut.problem; %0 is converged
            Qi{ii}=value(Q_sigI); gamma_i(ii)=value(gamma);

        end %end lmi

        %save topology optimization outputs
        convergenceTop{ss}=convergenceFlags;
        Qi_top{ss}=Qi; gamma_i_top{ss}=gamma_i;

    end %end topology

    %get all convergence flags
    allTopConverge=[convergenceTop{1}; convergenceTop{2}; convergenceTop{3}; convergenceTop{4}];
    
    %save off outputs for this tau
    allTauConvergence{iTau}=allTopConverge;

end %end tau sweep

%make convergence plots
top1=1:2; top2=3; top3=4; top4=5:6;
convergenceFlag=zeros(nTop,lenTau);
%parse convergence
for ttau=1:lenTau
    convergenceFlag(1,ttau)=sum(allTauConvergence{ttau}(top1));
    convergenceFlag(2,ttau)=sum(allTauConvergence{ttau}(top2));
    convergenceFlag(3,ttau)=sum(allTauConvergence{ttau}(top3));
    convergenceFlag(4,ttau)=sum(allTauConvergence{ttau}(top4));
end
[XX,YY]=meshgrid((1:nTop),tauVec);
figure('Name', 'Convergence');
surf(XX,YY,convergenceFlag'); colormap winter; view(2);
ylabel("\tau [sec]"); xlabel("topology"); colorbar;

%function to get one of four graph topologies as described in the paper
function [L_thetaBar,L_theta,d_sigI,nComponents]=GetGraphTopology(aij,iTop)
%get the graph matricies for one of the 4 topologies
%this is hardcoded to work with 6 agents for the examples, 
%dont pass in number of agents as an argument
nAgents=6; d_sigI=[];
[A_theta,E_sigma]=deal(zeros(nAgents));
%rows=talks to, cols=receives from
switch iTop
    case 1%a
        A_theta(1,6)=aij; A_theta(6,1)=aij;
        A_theta(3,4)=aij; A_theta(4,3)=aij;
        E_sigma(1,1)=1; E_sigma(6,2)=1;
        E_sigma(3,3)=1; E_sigma(4,4)=1;
        d_sigI(1)=2; d_sigI(2)=2;
    case 2%b
        A_theta(1,2)=aij; A_theta(2,1)=aij;
        E_sigma(1,1)=1; E_sigma(2,2)=1;
        d_sigI(1)=2;
    case 3%c
        A_theta(6,5)=aij; A_theta(5,6)=aij;
        A_theta(4,5)=aij; A_theta(5,4)=aij;
        E_sigma(4,1)=1; E_sigma(5,2)=1; E_sigma(6,3)=1;
        d_sigI(1)=3;
    case 4%d
        A_theta(2,3)=aij; A_theta(3,2)=aij;
        A_theta(5,6)=aij; A_theta(6,5)=aij;
        E_sigma(2,1)=1; E_sigma(3,2)=1;
        E_sigma(5,3)=1; E_sigma(6,4)=1;
        d_sigI(1)=2; d_sigI(2)=2;
end
%calculate number of components = l_sig
nComponents=length(d_sigI);
%calculate D
D_theta=diag(sum(A_theta,2));
%calculate L, paper uses non-normalized version
L_theta=(D_theta - A_theta);
%calculate the permutated graph laplacian    
L_thetaBar=E_sigma'*L_theta*E_sigma; 
end