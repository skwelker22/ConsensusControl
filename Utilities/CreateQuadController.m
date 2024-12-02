function [CLss, Kquad,nStates,nControls,nObs]=...
                CreateQuadController(dT,nDims,truthDynamics,plotFlag)
%dynamics : taken from Modeling, design and control of a 6-DOF quadcopter
%fleet ... by Anshuman
g=9.81;
% g=1;

%set helper variables to develop state space model
nStates=2*nDims; nControls=nDims; nObs=nDims;

%state space
A=zeros(nStates); A(1:nDims,nDims+1:nStates)=eye(nDims);
B=zeros(nStates,nControls); 
C=zeros(nObs,nStates); C(1:nDims,1:nDims)=eye(nDims);
D=zeros(nObs,nControls);

% set B matrix
switch nDims
    case 2
        B(3,1)=-g; B(4,2)=g;
    case 3
        B(4,1)=-g; B(5,2)=g; B(6,3)=1.5015;
    otherwise
        error(['controller for this case not implemented, ',...
              'please enter 2 or 3 for nDims']);
end

%create state space model
Pquad=ss(A,B,C,D); 

%if truth dynamics enabled, just output constant velocity dynamics
%with gravity induced B matrix
if(truthDynamics)    
    CLss=Pquad; Kquad=zeros(nControls,nStates);
else
    %% position controller design
    % muK=900; QK=2e-1;
    QK=0.001;
    % Kquad=-lqr(A,B,C'*C+A'*C'*C*A/muK,QK);
    Kquad=-lqr(A,B,C'*C,QK);

    %closed loop sensistivity and comp sensitivity
    kr=-inv(C*inv(A+B*Kquad)*B);
    CLss=ss(A+B*Kquad, B*kr, C, zeros(1,1));

    %if make plots
    if(plotFlag)
    % figure('Name', 'T,S Bode');
    % bodemag(CLss);
    figure('Name','Closed Loop Comp Sens Step');
    h=stepplot(CLss);
    h.showCharacteristic('RiseTime');
    end
end

end