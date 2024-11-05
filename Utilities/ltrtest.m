%% consensus control flocking simulation
%startup & format
clear all; close all; home;

format long g;

%format
set(0, 'DefaultTextFontName', 'Arial', 'DefaultTextFontSize', 20, 'DefaultTextFontWEight', 'Bold');
set(0, 'DefaultAxesFontName', 'Arial', 'DefaultAxesFontSize', 20, 'DefaultAxesFontWeight', 'Bold', 'DefaultAxesLineWidth', 1.5);
set(0, 'DefaultLineLineWidth', 3, 'DefaultLineMarkerSize', 10);

addpath(genpath('/Users/swelker/Desktop/swelker/ConsensusControl'));

%%
%a=[-10 1;0 -1], b=[0;1],c=[20 -1 ],d=0;
%P=ss(a,b,c,d);
%P=tf([-.1 1],[.1 1])*tf(1,[1 0 -1]);P=ss(P); 
%P=tf([1 -2],[1 -4 3]);P=ss(P);
%P=tf(1,[1 0 0 ]);
%P=tf(.1,conv([1 .1],[1 .1]));
%P=ss(P);
%quad dynamics
g=9.81;
nStates=2; nControls=1; nObs=1;
A=zeros(nStates);
A(1,2)=eye(nStates/2);
B=zeros(nStates,nControls);
B(2,1)=-g; %B(5,2)=g; B(6,3)=1.5015;
C=zeros(nObs,nStates);
C(1,1)=eye(nStates/2);
D=zeros(nObs,nControls);
P=ss(A,B,C,D);
%transfer function
%Ps=Cp*(s*eye(nStates)-Ap)^(-1)*Bp + Dp;

wn=10;wm=12;D1=tf([1 2*.1*wn wn^2],[1 2*.2*wm wm^2])*wm^2/wn^2;
P=P*ss(D1);
D2=ss(tf([-.05 1],[.05 1]));
D=D1*D2;
%P=P*D2
BW=5;

pid=pidqtune(BW/1.7,P,[],[1 1e-7 0]);C0=tf(pid(1,:),pid(2,:))

[pid2,C1]=pidpmtune(BW*1.5,P,.01,50)
%clf, 
figure,
bode(P*C0,P*C1)

S=fbk(1,P*C0);T=fbk(P*C0,1);

%% LQG gains

In = ss(tf(1,[1  0]));   Pa=P*In;  % int cascade at input
a=Pa.a;b=Pa.b;c=Pa.c;d=Pa.d;

L=-(lqr(a',c',1*eye(size(a))+1e12*b*b',1))'
K=-(lqr(a,b,c'*c+a'*c'*c*a/900,2.e-1))
%clf
S0=fbk(1,P*C0);
SS3=ss(a+b*K,b,K,1);
figure,bodemag(SS3,S0)
%ppc gains 1
% Qc=ctrb(a,b); K=-([0 1]*inv(Qc)*(a+10*eye(2,2))*(a+10*eye(2,2)));
% Qo=ctrb(a',c'); L=-([0 1]*inv(Qo)*(a'+10*eye(2,2))*(a'+10*eye(2,2)))';
%ppc gains 2
% Qc=ctrb(a,b); K=-([0 0 1]*inv(Qc)*(a+BW*eye(size(a)))*(a+BW*eye(size(a)))*(a+BW*eye(size(a))));
% Qo=ctrb(a',c'); L=-([0 0 1]*inv(Qo)*(a'+BW*eye(size(a)))*(a'+BW*eye(size(a)))*(a'+BW*eye(size(a))))';

%look at Glover/McFarlen controller


%%
C3=In*ss(a+b*K+L*c,L,K,0); tf(C3)
S3=fbk(1,P*C3);
figure,bodemag(S3,C3,S0,C0), %pause
figure,step(fbk(P*C3,1),fbk(P*C0,1))

% 
% L1=-(lqr(a',c',eye(size(a))+1e8*b*b',1))'
% K=-(lqr(a,b,c'*c+a'*c'*c*a/900,1e-2))
% C31=In*ss(a+b*K+L1*c,L1,K,0);tf(C31)
% S31=fbk(1,P*C31);
% SS3=ss(a+b*K,b,K,1);
% bodemag(S3,S31,SS3,C3,C31), title('integrator cascade at input recovery')
% pause
% step(fbk(P*C3,1),fbk(P*C31,1))


%  =IF(C3>=1000,"a+",IF(C3>=900,"A",IF(C3>=800,"b",IF(C3>=600,"C",IF(C3>=500,"D","?")))))

%  =IF(C3>=850,"a",IF(C3>=750,"B",IF(C3>=600,"C",IF(C3>=500,"D",IF(C3>=400, "DD","?")))))
% =IF(C3>=900,"a",IF(C3>=800,"B",IF(C3>=650,"C",IF(C3>=500,"D",IF(C3>=400,"DD","?")))))

% =IF(C3>=180,"A",IF(C3>=170,"a-",IF(C3>=160,"b+",IF(C3>=150,"B",IF(C3>=140,"b-","?")))))

%%
T=0.1;  
z=tf([1 0],[1],T) ; % define z as a transfer function
z=tf('z',T)       ; % MATLAB also understands this definition
s=(z-1)/T         ; % define the FE rule
HD = 1/(s^2-s+2)        % Write the transfer function in terms of s 
                   % The result in H is the DT FE equivalent
                   
T=0.1;  
s=tf([1 0],[1])   ; % define s as a transfer function
s=tf('s')         ; % MATLAB also understands this definition
z=1+s*T           ; % define the FE rule
HC = 1/(z^2-z+2)      % Write the transfer function in terms of z 
                   % The result in H is the CT FE equivalent

                   
T=0.001;  
s=tf([1 0],[1])   ; % define s as a transfer function
z=1+s*T           ; % define the FE rule
HC2 = .1/(z-.95)      % Write the transfer function in terms of z 
                   % The result in H is the CT FE equivalent
