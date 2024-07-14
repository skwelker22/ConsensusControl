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
nSec      = 50;
% dT        = 0.4;
dT        = 0.05;
nSamps    = floor(nSec/dT);
tt        = (0:dT:nSec-dT)';
makeMovie = 0;

%define graph dynamics
A = [0, 1; 0, 0];
B = [0, 1]';
F = [1, 0];

%design control value k
kGain = -1.6;
K     = [0, kGain];

%get graph parameters
nNodes = 10; %M from paper

%define switching network as random switch
graphIx = randi([1,4], nSamps, 1);
A1 = zeros(nNodes,nNodes);
A1(1,2) = 1;  A1(1,9)  = 1;  A1(1,10) = 1;
A1(2,1) = 1;  A1(2,3)  = 1;  A1(2,8)  = 1;
A1(3,2) = 1;  A1(3,4)  = 1;  A1(3,7)  = 1;
A1(4,3) = 1;  A1(4,5)  = 1;  A1(4,6)  = 1;
A1(5,4) = 1;  A1(5,6)  = 1;
A1(6,4) = 1;  A1(6,5)  = 1;  A1(6,7)  = 1;
A1(7,3) = 1;  A1(7,6)  = 1;  A1(7,8)  = 1;
A1(8,2) = 1;  A1(8,7)  = 1;  A1(8,9)  = 1;
A1(9,1) = 1;  A1(9,8)  = 1;  A1(9,10) = 1;
A1(10,1) = 1; A1(10,9) = 1;  

A2 = A1;
A2(1,9) = 0; A2(9,1) = 0;
A2(2,8) = 0; A2(8,2) = 0;
A2(3,7) = 0; A2(7,3) = 0;
A2(4,6) = 0; A2(6,4) = 0;

A3 = A1;
A3(2,8) = 0; A3(8,2) = 0;
A3(3,7) = 0; A3(7,3) = 0;

A4 = A1;
A4(1,9) = 0; A4(9,1) = 0;
A4(4,6) = 0; A4(6,4) = 0;

ANetworkType = cell(4,1);
ANetworkType{1} = A1; ANetworkType{2} = A2; 
ANetworkType{3} = A3; ANetworkType{4} = A4;

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

%loop over graphs
for nn = 1:(nSamps-1)

    %get current network configuration
    AGraph = ANetworkType{graphIx(nn)};

    %in degree
    DGraph = diag(sum(AGraph, 2));

    %graph laplacian
    LGraph = (DGraph - AGraph);
    LGraph_norm = DGraph^(-1)*LGraph;
    
    if makeMovie
        maxDELTA = 3;
        [eVecL, eValL] = eig(LGraph);
        DELTA = max(diag(DGraph));
        th = 0:pi/50:2*pi;
        gershGorinDiskX = DELTA * cos(th) + DELTA;
        gershGorinDiskY = DELTA * sin(th) + 0;
        figure(111);
        scatter(real(eValL), imag(eValL), 'k', 'LineWidth', 4);
        xlabel("Real{\lambda_L}"); ylabel("Imag{\lambda_L}");
        hold on; grid on;
        plot(gershGorinDiskX, gershGorinDiskY, '--b');
        ylim([-maxDELTA-1,maxDELTA+1]); xlim([-1, 2*maxDELTA+1])
        drawnow;
        pause(0.1);
        hold off;
    end

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

%% plots
ccc = turbo(nNodes);
figure('Name', 'Fixed Topology');
colIx = 1; 
posIx = [1:2:2*nNodes];
velIx = [2:2:2*nNodes];
for ii = 1:nNodes
    subplot(121);
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

%graph plot
figure('Name', 'Network Configuration');
plot( tt, graphIx, '-ob' );
xlabel("Time [sec]"); ylabel("Topology");
yticks([1:4]);
yticklabels({'G_1', 'G_2', 'G_3', 'G_4'});
