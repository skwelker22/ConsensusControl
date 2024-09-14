%start
clear all; close all; home;

%format
set(0, 'DefaultTextFontName', 'Arial', 'DefaultTextFontSize', 20, 'DefaultTextFontWEight', 'Bold');
set(0, 'DefaultAxesFontName', 'Arial', 'DefaultAxesFontSize', 20, 'DefaultAxesFontWeight', 'Bold', 'DefaultAxesLineWidth', 1.5);
set(0, 'DefaultLineLineWidth', 3, 'DefaultLineMarkerSize', 10);
format long g;

%create quad object
m = 0.1;         %kg
Ixx = 0.00062;   %kg-m^2
Iyy = 0.00113;   %kg-m^2
Izz = 0.9*(Ixx + Iyy); %kg-m^2 (Assume nearly flat object, z=0)
dx = 0.114;      %m
dy = 0.0825;     %m

%configuration structure of some physics related parameters for the quad
quadConfig=struct('m', m, 'Ixx', Ixx, 'Iyy', Iyy, 'Izz', Izz, ...
                  'dx', dx, 'dy', dy);

% Simulation time and model parameters
tstep = 0.02;          % Sampling time (sec)
simulation_time = 30;   % Length of time to run simulation (sec)
t = (0:tstep:simulation_time);   % time array
nTime=length(t);

% Model size
n_states = 12; n_inputs = 4;

% Initialize State Conditions
x = zeros(n_states,nTime);
%initial height
x(12,1) = 0.0;

% Initialize inputs
u = zeros(n_inputs,nTime);  % time history of input vectors
% Initial control inputs
u(:,1) = zeros(4,1);

%create quad object
quad=Quadcopter(quadConfig);

% March through time array and numerically solve for vehicle states
nSecDebug=5;
for k = 1:(nTime-1)
    %display time every second
    if(mod(k,nSecDebug/tstep)==0)
        fprintf('t= %3.1f [sec]\n', k*tstep);
    end

    % Determine control inputs based on current state
    u(:,k) = controlInputs(x(:,k), t(k));
    
    %propogate dynamics
    x(:,k+1) = quad.PropogateDynamics(x(:,k), u(:,k), tstep);
end

figure('Name','States');
subplot(311);
plot(t,x(12,:),'b');
ylabel('h (m)'); title('Time History of Height, X Position, and Pitch');
subplot(312);
plot(t,x(10,:),'b');
ylabel('x (m)');
subplot(313);
plot(t,x(8,:)*Quadcopter.RTD,'b');
ylabel('Theta (deg)'); xlabel('Time (sec)');

figure('Name', 'Horiz distance vs. Height');
plot(x(10,1:20:end),x(12,1:20:end),'bo-'); hold on;
plot(x(10,1) + 0.1, x(12,1),'x');
plot(x(10,end), x(12,end),'o');
ylabel('h [m]'); xlabel('x [m]');
title('Vertical Profile'); axis equal;

figure('Name', 'Control Input');
plot(t,u(1,:),'b'); hold on; grid on;
plot(t,u(2,:),'g');
plot(t,u(3,:),'r');
plot(t,u(4,:),'k');
xlabel('Time (sec)'); ylabel('Propeller RPM');
title('Time History of Control Inputs');
legend('T1', 'T2', 'T3', 'T4');

function u=controlInputs(x, t)
% Inputs: Current state x[k], time t
% Returns: Control inputs u[k]
% Placeholder Function ####

% Trim RPM for all 4 propellers to provide thrust for a level hover
trim = 3200;

pitch_cmd = 0;
roll_cmd = 0;
climb_cmd = 0;
yaw_cmd = 0;

% Example open loop control inputs to test dynamics:
%  Climb
if t < 11.0
    climb_cmd = 500;
end

%  Pitch Forward
if t > 8.0
    pitch_cmd = -10;
end
if t > 9.0
    pitch_cmd = 10;
end
if t > 10.0
    pitch_cmd = 0;
end

%  Pitch Backward
if t > 12.0
    pitch_cmd = 15;
end
if t > 13.0
    pitch_cmd = -15;
end
if t > 14.0
    pitch_cmd = 0;
end

% Increase lift
if t > 16.0
    climb_cmd = 150;
end

% RPM command based on pitch, roll, climb, yaw commands
u = zeros(4,1);
u(1) = trim + ( pitch_cmd + roll_cmd + climb_cmd - yaw_cmd) / 4;
u(2) = trim + (-pitch_cmd - roll_cmd + climb_cmd - yaw_cmd) / 4;
u(3) = trim + ( pitch_cmd - roll_cmd + climb_cmd + yaw_cmd) / 4;
u(4) = trim + (-pitch_cmd + roll_cmd + climb_cmd + yaw_cmd) / 4;

end


