%% consensus control flocking simulation
%startup & format
clear all; close all; home;

format long g;

%format
set(0, 'DefaultTextFontName', 'Arial', 'DefaultTextFontSize', 20, 'DefaultTextFontWEight', 'Bold');
set(0, 'DefaultAxesFontName', 'Arial', 'DefaultAxesFontSize', 20, 'DefaultAxesFontWeight', 'Bold', 'DefaultAxesLineWidth', 1.5);
set(0, 'DefaultLineLineWidth', 3, 'DefaultLineMarkerSize', 10);

%add lmi paths
addpath(genpath('C:\Users\skwel\sam\School\MAE509_LMIOPT\YALMIP-master'));
addpath('C:\Program Files\Mosek\10.2\toolbox\r2017a');

%% setup plant and design problem
nStates=7; %[x, y, z, x_dot, y_dot, z_dot, psi]
nControls=4; %[theta, phi, T, r]
nMeas=4; %[x,y,z,psi]
nAugStates=4;
g=9.81;

%define the plant dynamics
Ap=zeros(nStates); 
Ap(1:3,4:6)=eye(3);

Bp=zeros(nStates,nControls); 
Bp(4,1)=-g; Bp(5,2)=g; Bp(6,3)=1.5015; Bp(7,4)=1;

Cp=zeros(nMeas,nStates);
Cp(1:3,1:3)=eye(3); Cp(4,7)=1;

s=tf('s');
Ps=Cp*(s*eye(nStates)-Ap)^(-1)*Bp;

%zero out near zero terms
Ps(1,1).Denominator{1}=round(Ps(1,1).Denominator{1});
Ps(2,2).Denominator{1}=round(Ps(2,2).Denominator{1});
Ps(3,3).Denominator{1}=round(Ps(3,3).Denominator{1});

%design servo controller

Ag=[Ap zeros(nStates,nControls); Cp, zeros(nControls,nControls)];
Bg=[Bp; zeros(nControls,nControls)];
Cg=eye(nStates);

%lq tuning
rho=0.2; R=rho.*eye(nControls,nControls);
Q=diag([10,10,10,10,10,10,100,100,100,100,1]);
N=zeros(nStates+nControls,nControls);

%run lqr to solve the riccati
[Kopt,Sopt]=lqr(Ag,Bg,Q,R,N);
Ki=Kopt(:,1:nAugStates); Ky=Kopt(:,nAugStates+1:end); Cr=eye(nAugStates);

%form the closed loop response
Acl=[Ap-Bp*Ky*Cp Bp*Ki; -Cr*Cp, zeros(nAugStates,nAugStates)];

