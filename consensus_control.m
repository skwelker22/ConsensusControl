%Consensus Control for a Class of Networks of Dynamic Agents
%% hw 3
%startup & format
clear all; close all; home;

format long g;

%format
set(0, 'DefaultTextFontName', 'Arial', 'DefaultTextFontSize', 20, 'DefaultTextFontWEight', 'Bold');
set(0, 'DefaultAxesFontName', 'Arial', 'DefaultAxesFontSize', 20, 'DefaultAxesFontWeight', 'Bold', 'DefaultAxesLineWidth', 1.5);
set(0, 'DefaultLineLineWidth', 3, 'DefaultLineMarkerSize', 10);

%% start algorithm
%iterations
nSec   = 20;
dT     = 0.4;
nSamps = floor(nSec/dT);
tt     = (0:dT:nSec-dT)';

%define graph dynamics
A = [0, 1; 0, 0];
B = [0, 1]';
F = [1, 0];

%design control value k
kGain = -1.6;
K     = [0, kGain];

% %define for graph G1
AGraph = [0, 1, 0, 1, 0, 0;
          1, 0, 1, 0, 0, 0;
          0, 1, 0, 0, 1, 0;
          1, 0, 0, 0, 1, 1;
          0, 0, 1, 1, 0, 1;
          0, 0, 0, 1, 1, 0];

%in degree
DGraph = diag(sum(AGraph, 2));

%graph laplacian
LGraph = (DGraph - AGraph);
LGraph_norm = DGraph^(-1)*LGraph;

%get graph parameters
nNodes = size(AGraph, 1); %M from paper
nEdges = 7; %counting

%create phi
IM  = eye(nNodes);
PHI = kron(IM, (A + B*K)) - kron(LGraph_norm, (B*F));

%get transfer function for the plant
Ps = ss(A, B, F, 0);
[Pb, Pa] = ss2tf(Ps.a, Ps.b, Ps.c, Ps.d);
PsDoubleIntegrator = tf(Pb, Pa);

%get closed loop with state feedback
nStates = 2; nControls = 1;
Ls_ss = ss(A+B*K, B, F, 0);
[Ls_b, Ls_a] = ss2tf(Ls_ss.a, Ls_ss.b, Ls_ss.c, Ls_ss.d);
Ls = tf(Ls_b, Ls_a);

%with PM = 45 degrees, wc = 1 rad/sec, T = 0.1 seconds,
T  = 0.0; omegaC = 5; PM = 45;
Ps = exp(-tf([1,0],1)*T)*PsDoubleIntegrator;
% C_PD = DesignPdController(Ps, omegaC, PM);
% Ls  =  Ps * C_PD;
[reLs, imLs, omegaLs] = nyquist(Ls);
reLs = squeeze(reLs(1,:,:)); imLs = squeeze(imLs(1,:,:));
th = 0:pi/50:2*pi;
unitCircX = cos(th);
unitCircY = sin(th);

%check spectral convergence conditions
[eVecL, eValL] = eig(LGraph);
eValL = diag(eValL);
DELTA = max(diag(DGraph));
gershGorinDiskX = DELTA * cos(th) + DELTA;
gershGorinDiskY = DELTA * sin(th) + 0;
smallEig = 1e-10;
nonZeroEigs = eValL(abs(eValL) > smallEig);
lenNonZeroEigs = length(nonZeroEigs);
zeroImagEigsIx = abs(imag(nonZeroEigs))<1e-16;
if (1)
    figure('Name', 'Nyquist with Laplacian');
    plot( reLs, imLs, 'k', 'LineWidth', 1 ); hold on; grid on;
    plot( reLs, -imLs, 'k', 'LineWidth', 1);
    plot( unitCircX, unitCircY, '--b', 'LineWidth', 1);
    if (sum(zeroImagEigsIx) == lenNonZeroEigs)
        plot(-1./real(nonZeroEigs), zeros(lenNonZeroEigs,1), 'xk', 'LineWidth', 4);
    else
        %this plotting is only correct when ALL imaginary parts of eigen
        %values are zero, otherwise need to do something different here.
        plot(-1./real(nonZeroEigs), -1./imag(nonZeroEigs), 'xk', 'LineWidth', 4);
    end
    xlim([-2,2]); ylim([-2,2]);
    xlabel("Real Axis"); ylabel("Imaginary Axis");

    figure('Name', 'Eigen Values of Laplacian');
    scatter(real(eValL), imag(eValL), 'k', 'LineWidth', 4);
    xlabel("Real{\lambda_L}"); ylabel("Imag{\lambda_L}");
    hold on; grid on;
    plot(gershGorinDiskX, gershGorinDiskY, '--b');
end

%form C (incidence matrix) NOTE: this C matrix only applies to the defined
%A
CGraph = zeros(size(LGraph));
CGraph(1,1) = -1; CGraph(1,4) = -1;
CGraph(2,1) =  1; CGraph(2,2) = -1;
CGraph(3,2) =  1; CGraph(3,3) = -1;
CGraph(4,4) =  1; CGraph(4,5) = -1; CGraph(4,6) = -1;
CGraph(5,3) =  1; CGraph(5,5) =  1; CGraph(5,7) = -1;
CGraph(6,6) =  1; CGraph(6,7) =  1;
incidenceEqualLaplacian = CGraph*CGraph' == LGraph;
incidenceFlag1 = sum(incidenceEqualLaplacian(:)) == numel(LGraph);
incidenceFlag2 = rank(LGraph) == rank(CGraph);

%both of the above properties need to be true for the graph
if ~(incidenceFlag1 && incidenceFlag2)
    warning('C*C^T != L, check your graph');
end

%check if graph is connected
isConnected = rank(LGraph) == (nNodes-1);
if isConnected
    disp("Graph is connected, rank(L) = M-1");
else
    warning("Graph is not fully connected");
end

%check eigen values of phi and property lim t->inf exp(phi * t) = wrwl^T
[eVecR, eVal, eVecL] = eig(PHI); %eig returns the right eigenvectors
eVal = diag(eVal);
wr = eVecR(:,abs(eVal) <1e-13);
wl = eVecL(:,abs(eVal) <1e-13);
wrTrue = 1/sqrt(nNodes) * kron(ones(nNodes,1), [1,0]');
wlTrue = 1/(-K(2)*sqrt(nNodes)) * kron(ones(nNodes,1), [-K(2), 1]');
wrTranwlTheory = wrTrue'*wlTrue;
wrTranwl = wr' * wl;

% if wrTranwl ~= 1
%     error('Left/Right eigen vector check failed');
% end
rlEigMatTrue = wrTrue*wlTrue';
rlEigMatEst = zeros(2*nNodes, 2*nNodes, nSamps);
for ii = 1:nSamps
    rlEigMatEst(:,:,ii) = expm(PHI * ii * dT);
end

if(0)
    figure('Name', 'Connected Graph PHI Check');
    plot( diag(rlEigMatTrue) );
    hold on; grid on;
    plot( diag(rlEigMatEst(:,:,end))./nNodes, '--' );
end


%initialize node states]
iLow  = -100; iHigh = 100;
x0    = randi([iLow, iHigh], nNodes, 1);
x0Avg = mean(x0);
v0    = zeros(nNodes,1);

%stack states
xi0 = [x0, v0]';
xi0 = xi0(:);

%propagate the consensus control law at each node
nXi        = length(xi0); 
xi         = zeros(nXi, nSamps);
nConsensus = nNodes - 1;
xi(:,1)    = xi0;
nodeIx     = reshape((1:2*nNodes)', 2, []);

runSerial = true;
if (runSerial)
    for nn = 1:(nSamps-1)
        %Laplacian version
        for ii = 1:nNodes
            %get self states
            nodeii = nodeIx(:,ii);
            xi_nn_self = xi(nodeii, nn);

            %form the consensus term
            xi_consensus = 0;
            for jj = 1:nNodes
                nodejj = nodeIx(:,jj);
                xi_nn_other = xi(nodejj, nn);
                xi_consensus = xi_consensus + AGraph(ii,jj) * (xi_nn_other - xi_nn_self);
            end

            %update states
            % xi(nodeii, nn+1) = xi_nn_self + xi_consensus./nConsensus; %works
            xi(nodeii, nn+1) = (eye(size(A)) + A + B*K) * xi_nn_self + B * F * xi_consensus./nConsensus;
        end
    end
else % run fsolve
    %continuous time version
    options = optimoptions('fsolve','Display','off', 'TolFun', 1E-8, 'TolX', 1E-8);
    for nn = 1:(nSamps-1)
        xi(:,nn+1) = fsolve(@(x)(eye(size(PHI)) + PHI)*x, xi(:,nn), options);
        % xi(:,nn+1) = fsolve(@(x)PHI*x, xi(:,nn) + sqrt(0.001)*randn(2*nNodes,1), options);
    end
end

%% plots
ccc = turbo(nNodes);
figure('Name', 'Fixed Topology');
colIx = 1; 
posIx = [1:2:2*nNodes];
velIx = [2:2:2*nNodes];
for ii = 1:nNodes
    subplot(121);
    % plot( tt, xi(posIx(ii),:), 'Color', ccc(colIx,:) );
    plot( tt, xi(posIx(ii),:), 'b', 'LineWidth', 3 );
    hold on; %grid on;
    xlabel("Time [sec]"); ylabel("Position");

    subplot(122);
    plot( tt, xi(velIx(ii),:), 'r', 'LineWidth', 3 );
    hold on; %grid on;
    xlabel("Time [sec]"); ylabel("Velocity");
    colIx = colIx + 1;
end
%plot average
subplot(121); yline(x0Avg, 'LineWidth', 3);

