%% consensus control flocking simulation
%startup & format
clear all; close all; home;

format long g;

%format
set(0, 'DefaultTextFontName', 'Arial', 'DefaultTextFontSize', 20, 'DefaultTextFontWEight', 'Bold');
set(0, 'DefaultAxesFontName', 'Arial', 'DefaultAxesFontSize', 20, 'DefaultAxesFontWeight', 'Bold', 'DefaultAxesLineWidth', 1.5);
set(0, 'DefaultLineLineWidth', 3, 'DefaultLineMarkerSize', 10);

%add libraries for controller design
if ispc
    addpath(genpath('C:\Users\skwel\sam\Research\ConsensusControl'));
else
    addpath(genpath('/Users/swelker/Desktop/swelker/ConsensusControl'));
end

%% test dcm
Rx=@(phi)[1,      0,        0;
          0, cosd(phi), sind(phi);
          0, sind(phi), cosd(phi)];
Ry=@(theta)[cosd(theta), 0, sind(theta);
                 0       1,      0     ;
           -sind(theta), 0, cosd(theta)];
Rz=@(xi)[cosd(xi), sind(xi), 0;
        -sind(xi), cosd(xi), 0;
            0         0      1];
Rzyx=@(phi,theta,xi)Rz(xi)*Ry(theta)*Rx(phi);

ux=[1,0,0]'; uy=[0,1,0]'; uz=[0,0,1]';

%% driver for controller
%controller settings
plotFlag=1; truthDynamics=0;

%get simulation configuration
rndSeed=3; ccfg = GetConsensusConfig(rndSeed);

%create controller
[AgentSS, ~]=CreateQuadController(ccfg.dT,ccfg.nDims,truthDynamics,plotFlag);

%set initial positions and velocities
qi0  = sqrt(ccfg.distVar)*randn(ccfg.nDims,ccfg.nNodes); %initial agent positions
[piHdgi,piHdgi2] = deal(0); %heading angle
pRngX = [-2,2]; pRngY = [-1,1]; pRngZ = [-1,1];
pi0X = randi(pRngX, 1, ccfg.nNodes);
pi0Y = randi(pRngY, 1, ccfg.nNodes);
pi0Z = randi(pRngZ, 1, ccfg.nNodes);

switch ccfg.nDims
    case 2
        pi0  = repmat([cos(piHdgi); sin(piHdgi)],1,ccfg.nNodes);
        x_star=[8;4]; v_star=[0;0];
    case 3
        pi0  = [pi0X; pi0Y; pi0Z];
        x_star=[15;0;0]; v_star=[0;0;0];
end

%set control to 0
uAlpha=zeros(ccfg.nDims,ccfg.nNodes,ccfg.nSamps);

%init member positions and velocities over time
[xi, vi, xi2, vi2] = deal(zeros(ccfg.nDims,ccfg.nNodes,ccfg.nSamps));
xi(:,:,1) = qi0; vi(:,:,1) = pi0;
xi2(:,:,1)=qi0; vi2(:,:,1)=pi0;

%discrete time system
AgentSSd=c2d(AgentSS, ccfg.dT);

%prop nav gain
kPN=4;

%set initial conditions
[rmtilast,rmti2last]=deal(zeros(ccfg.nDims,ccfg.nNodes));
for ii = 1:ccfg.nNodes
rmtilast(:,ii)=x_star-xi(:,ii,1); 
rmti2last(:,ii)=x_star-xi2(:,ii,1);
end

%propogate ahead with constant velocity
xi(:,:,1)=xi(:,:,1)+vi(:,:,1)*ccfg.dT;
xi2(:,:,1)=xi2(:,:,1)+vi2(:,:,1)*ccfg.dT;

%propagate dynamics with dt and ct systems
simEnd=false;
for tt = 1:ccfg.nSamps
    for ii = 1:ccfg.nNodes
        %relative position
        rmti = x_star-xi(:,ii,tt); rmti2=x_star-xi2(:,ii,tt);
        
        %range to target
        Ri = norm(rmti); Ri2=norm(rmti2);

        %relative velocity
        vmti=(rmti-rmtilast(:,ii))./ccfg.dT; vmti2=(rmti2-rmti2last(:,ii))./ccfg.dT;

        %agent velocity
        vpi=norm(vi(:,ii,tt)); vpi2=norm(vi2(:,ii,tt));

        %los angle
        lambdai=atan2(rmti(2),rmti(1)); lambdai2=atan2(rmti2(2),rmti2(1));

        %los angle rate
        lambdaDoti=(rmti(1)*vmti(2) - rmti(2)*vmti(1)) / Ri^2;
        lambdaDoti2=(rmti2(1)*vmti2(2) - rmti2(2)*vmti2(1)) / Ri2^2;
        
        %calculate the command
        uki=kPN*vpi*lambdaDoti;
        uki2=kPN*vpi2*lambdaDoti2;

        % Terminate sim at intercept
        if (abs(Ri) <= 0.25 || abs(Ri2) <=0.25)
            ttFinal=tt; simEnd=true;
            break;
        end

        %calculate commands
        uki=[-uki*sin(piHdgi); uki*cos(piHdgi)]; uki2=[-uki2*sin(piHdgi2); uki2*cos(piHdgi2)];

        xk=[xi(:,ii,tt); vi(:,ii,tt)];
        xkplus=AgentSSd.A*xk + AgentSSd.B*uki;
        xi(:,ii,tt+1)=xkplus(1:ccfg.nDims); vi(:,ii,tt+1)=xkplus(ccfg.nDims+1:end);

        % constant velocity system
        xi2(:,ii,tt+1) = xi2(:,ii,tt) + vi2(:,ii,tt)*ccfg.dT;
        vi2(:,ii,tt+1) = vi2(:,ii,tt) + uki2*ccfg.dT;

        %update heading
        piHdgi=atan2(vi(2,ii,tt+1),vi(1,ii,tt+1));
        piHdgi2=atan2(vi2(2,ii,tt+1),vi2(1,ii,tt+1));
        
        %update last
        rmtilast(:,ii)=rmti; rmti2last(:,ii)=rmti2;
    end
    if simEnd, break, end
end

%plot trajectories
ccplt=summer(ccfg.nNodes); ccplt2=winter(ccfg.nNodes);
figure('Name', 'Position Trajectories'); 
for tt=1:(ttFinal-1)
for ii=1:ccfg.nNodes
if(tt==1)
plot(xi(1,ii,tt), xi(2,ii,tt), 'kx', 'linewidth', 5); hold on;
plot(xi2(1,ii,tt), xi2(2,ii,tt), 'kx', 'linewidth', 5);hold on;
else
plot(xi(1,ii,tt), xi(2,ii,tt), 'marker', 'o', 'Color', ccplt(ii,:)); hold on;
plot(xi2(1,ii,tt), xi2(2,ii,tt), 'marker', 'o', 'Color', ccplt2(ii,:), 'linewidth', 2);hold on;
end
end
end
xlabel("East [m]"); ylabel("North [m]");

function [vCross]=TwoDimCross(v1,v2)
    vCross=v1(1)*v2(2)-v1(2)*v2(1);
end
