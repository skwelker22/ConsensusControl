function [Tquad_ss, Kquad,nStates,nControls,nObs]=...
                CreateQuadController(dT,nDims,truthDynamics,plotFlag)
%dynamics : taken from Modeling, design and control of a 6-DOF quadcopter
%fleet ... by Anshuman
g=9.81;
%g=1;

%set helper variables to develop state space model
nStates=2*nDims; nControls=nDims; nObs=nDims;

%revisit these dynamics because they appear to be WRONG***
A=eye(nStates); A(1:nDims,nDims+1:nStates)=dT*eye(nDims);
B=zeros(nStates,nControls); 
C=zeros(nObs,nStates); C(1:nDims,1:nDims)=eye(nDims);
D=zeros(nObs,nControls);

% set B matrix
switch nDims
    case 2
        B(3,1)=-g; B(4,2)=g;
        %B(3,1)=g; B(4,2)=g;
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
    Tquad_ss=Pquad; Kquad=zeros(nControls,nStates);
else
    %% controller design
    muK=900; QK=2e-1;
    Kquad=-lqr(A,B,C'*C+A'*C'*C*A/muK,QK);

    %closed loop sensistivity and comp sensitivity
    Squad_ss=ss(A+B*Kquad,B,Kquad,ones(nObs,nControls)); Tquad_ss=1-Squad_ss;

    %if make plots
    if(plotFlag)
    figure('Name', 'T,S Bode');
    bodemag(Squad_ss,Tquad_ss);
    nSec=50; tVec=[0:dT:nSec-dT];
    figure('Name','Closed Loop Comp Sens Step');
    step(tVec,Tquad_ss);
    end
end

end