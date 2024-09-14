classdef Quadcopter < handle

    properties
        % Physical Constants
        % source: https://gist.github.com/charlestytler
        m;   %kg
        Ixx; %kg-m^2
        Iyy; %kg-m^2
        Izz; %kg-m^2 (Assume nearly flat object, z=0)
        dx;  %m
        dy;  %m
    end

    properties(Constant)
        g = 9.81  %m/s/s
        RTD = 57.3; 
        DTR = 1/Quadcopter.RTD;
        optimzerOptions = optimoptions('fsolve', 'Display', 'off');
    end
    
    methods
        %class constructor
        function obj=Quadcopter(config)
            obj.m   = config.m;   
            obj.Ixx = config.Ixx; 
            obj.Iyy = config.Iyy; 
            obj.Izz = config.Izz; 
            obj.dx  = config.dx;  
            obj.dy  = config.dy;  
        end

        %propogate dynamics
        function x_dot=PropogateDynamics(obj,x,u,dt)
            x_dot=obj.RK4(x,u,dt);
        end
    end

    %static methods
    %don't require an instance of the class to call
    methods(Access=private)
        function residual = thrustEqn(obj, vi, prop_params)
            % Unpack parameters
            R=prop_params.R; A=prop_params.A; rho=prop_params.rho;
            a=prop_params.a; b=prop_params.b; c=prop_params.c;
            eta=prop_params.eta; theta0=prop_params.theta0; theta1=prop_params.theta1;
            U=prop_params.U; V=prop_params.V; W=prop_params.W; Omega=prop_params.Omega;
            
            % Calculate local airflow velocity at propeller with vi, V'
            Vprime = sqrt(U^2 + V^2 + (W - vi)^2);
            
            % Calculate Thrust averaged over one revolution of propeller using vi
            Thrust = 1/4 * rho * a * b * c * R * ...
                ( (W - vi) * Omega * R + 2/3 * (Omega * R)^2 * (theta0 + 3/4 * theta1) + ...
                  (U^2 + V^2) * (theta0 + 1/2 * theta1) );
            
            % Calculate residual for equation: Thrust = mass flow rate * delta Velocity
            residual = eta * 2 * vi * rho * A * Vprime - Thrust;
        end
        
        function Thrust=Fthrust(obj, x, u, dx, dy)
            % Inputs: Current state x[k], Commanded Propeller RPM inputs u[k],
            %         Propeller location distances dx, dy (m)
            % Returns: Thrust vector for 4 propellers (Newtons)
            % Propeller Configuration parameters
            R = 0.0762;     % propeller length/ disk radius (m) 
            A = pi * R ^ 2;
            rho = 1.225;    %kg/m^3  at MSL
            a = 5.7;        % Lift curve slope used in example in Stevens & Lewis
            b = 2;          % number of blades
            c = 0.0274;     % mean chord length (m)
            eta = 1;        % propeller efficiency
            
            % Manufacturer propeller length x pitch specification:
            p_diameter = 6; %inches
            p_pitch = 3;    %inches
            
            theta0 = 2*atan2(p_pitch, (2 * pi * 3/4 * p_diameter/2));
            theta1 = -4 / 3 * atan2(p_pitch, 2 * pi * 3/4 * p_diameter/2);

            % Local velocity at propeller from vehicle state information
            ub=x(1); vb=x(2); wb = x(3);
            p=x(4); q=x(5); r = x(6);

            % Transform velocity to local propeller location:
            %     [U,V,W] = [ub,vb,wb] + [p,q,r] x [dx,dy,0]
            U = ub - r * dy;
            V = vb + r * dx;
            W = wb - q * dx + p * dy;
            
            % Convert commanded RPM to rad/s
            Omega = 2 * pi / 60 * u;
            
            %Collect propeller config, state, and input parameters
            prop_params = struct('R',R,'A',A,'rho',rho,'a',a,'b',b,...
                                 'c',c,'eta',eta,'theta0',theta0,'theta1',theta1,...
                                 'U',U,'V',V,'W',W,'Omega',Omega);

            %define cost for optimization
            optFunc=@(vi) obj.thrustEqn(vi,prop_params);
            
            % Numerically solve for propeller induced velocity, vi
            % using nonlinear root finder, fsolve, and prop_params
            vi0 = 0.1;    % initial guess for vi
            vi = fsolve(optFunc,vi0,Quadcopter.optimzerOptions);
            
            % Plug vi back into Thrust equation to solve for T
            Vprime = sqrt(U^2 + V^2 + (W - vi)^2);
            Thrust = eta * 2 * vi * rho * A * Vprime;
        end
        
        %quad dynamics
        function xdot=stateDerivative(obj,x,u)
            % Inputs: state vector (x), input vector (u)
            % Returns: time derivative of state vector (xdot)
            %  State Vector Reference:
            %idx  0, 1, 2, 3, 4, 5,  6,   7,   8,   9, 10, 11
            %x = [u, v, w, p, q, r, phi, the, psi, xE, yE, hE]
            
            % Store state variables in a readable format
            ub = x(1);
            vb = x(2);
            wb = x(3);
            p = x(4);
            q = x(5);
            r = x(6);
            phi = x(7);
            theta = x(8);
            psi = x(9);
            xE = x(10);
            yE = x(11);
            hE = x(12);
            
            % Calculate forces from propeller inputs (u)
            F1 = obj.Fthrust(x, u(1),  obj.dx,  obj.dy);
            F2 = obj.Fthrust(x, u(2), -obj.dx, -obj.dy);
            F3 = obj.Fthrust(x, u(3),  obj.dx, -obj.dy);
            F4 = obj.Fthrust(x, u(4), -obj.dx,  obj.dy);

            Fz = F1 + F2 + F3 + F4;
            L = (F1 + F4) * obj.dy - (F2 + F3) * obj.dy;
            M = (F1 + F3) * obj.dx - (F2 + F4) * obj.dx;
            % N = -T(F1,obj.dx,obj.dy) - T(F2,obj.dx,obj.dy) + ...
            %     T(F3,obj.dx,obj.dy) + T(F4,obj.dx,obj.dy);
            %placeholder for now, need to implement a method that returns torque about the CG
            N = 0; 
            
            % Pre-calculate trig values
            cphi = cos(phi);   sphi = sin(phi);
            cthe = cos(theta); sthe = sin(theta);
            cpsi = cos(psi);   spsi = sin(psi);
            
            % Calculate the derivative of the state matrix using EOM
            xdot = zeros(12,1);
            
            xdot(1) = -Quadcopter.g * sthe + r * vb - q * wb;  % = udot
            xdot(2) = Quadcopter.g * sphi*cthe - r * ub + p * wb; % = vdot
            xdot(3) = 1/obj.m * (-Fz) + Quadcopter.g*cphi*cthe + q * ub - p * vb; % = wdot
            xdot(4) = 1/obj.Ixx * (L + (obj.Iyy - obj.Izz) * q * r);  % = pdot
            xdot(5) = 1/obj.Iyy * (M + (obj.Izz - obj.Ixx) * p * r);  % = qdot
            xdot(6) = 1/obj.Izz * (N + (obj.Ixx - obj.Iyy) * p * q);  % = rdot
            xdot(7) = p + (q*sphi + r*cphi) * sthe / cthe;  % = phidot
            xdot(8) = q * cphi - r * sphi;  % = thetadot
            xdot(9) = (q * sphi + r * cphi) / cthe;  % = psidot
            
            xdot(10) = cthe*cpsi*ub + (-cphi*spsi + sphi*sthe*cpsi) * vb + ...
                (sphi*spsi+cphi*sthe*cpsi) * wb;  % = xEdot
            
            xdot(11) = cthe*spsi * ub + (cphi*cpsi+sphi*sthe*spsi) * vb + ...
                (-sphi*cpsi+cphi*sthe*spsi) * wb; % = yEdot
            
            xdot(12) = -1*(-sthe * ub + sphi*cthe * vb + cphi*cthe * wb); % = hEdot
        end

        % 4th Order Runge Kutta Calculation
        function x_next=RK4(obj,x,u,dt)
            % Inputs: x[k], u[k], dt (time step, seconds)
            % Returns: x[k+1]

            % Calculate slope estimates
            K1 = obj.stateDerivative(x, u);
            K2 = obj.stateDerivative(x + K1 * dt / 2, u);
            K3 = obj.stateDerivative(x + K2 * dt / 2, u);
            K4 = obj.stateDerivative(x + K3 * dt, u);

            % Calculate x[k+1] estimate using combination of slope estimates
            x_next = x + 1/6 * (K1 + 2*K2 + 2*K3 + K4) * dt;
        end

    end
end