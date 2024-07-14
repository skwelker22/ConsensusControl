%% consensus control flocking simulation
%startup & format
clear all; close all; home;

format long g;

%format
set(0, 'DefaultTextFontName', 'Arial', 'DefaultTextFontSize', 20, 'DefaultTextFontWEight', 'Bold');
set(0, 'DefaultAxesFontName', 'Arial', 'DefaultAxesFontSize', 20, 'DefaultAxesFontWeight', 'Bold', 'DefaultAxesLineWidth', 1.5);
set(0, 'DefaultLineLineWidth', 3, 'DefaultLineMarkerSize', 10);

%% sim parms
%sim update rate
dT = 0.01;
dF = 1/dT;
tFinal = 100; %seconds
nSamps = tFinal/dT;
d = 7; %constraint distance agent to agent
r = 1.2 * d; %ball radius which defines the max distance for the proximity net
d_prime = 0.6 * d; %constraint distance between agents and obstacles
r_prime = 1.2 * d_prime; %ball radius defining an obstacle as a neigbor to an agent
eps_parm = 0.1; %fixed epsilon used for sigma-norm calculations
a = 5; %a parameter for potential function
b = a; %b parameter for potential function
c = abs(a-b)/sqrt(4*a*b); % c parameter for potential function
c1 = 1; c2 = 1; %c1,c2 for navigational feedback control term
% distVar = 2500; %variance on the intial positions of the nodes
distVar = 500;
nNodes  = 50; %number of  agents in the networks
nDims   = 3; %number of dimensions of the problem
qi0 = sqrt(distVar)*randn(nDims,nNodes); %initial agent positions
pRng = [-2,-1].^2; %range of initial agent velocities
pi0 = randi(fliplr(pRng), nDims, nNodes); %initial velocity
pi_angle = atan2d(pi0(2,:), pi0(1,:)); %relative heading angle

%define static gamma agent
qd = [[200; 30; 100], [50;40;10]]; pd = [5;0;0];

%define switch times
nSwitch = 10; switchIx = 1;
tNavSwitch = unique(randi([10,tFinal], nSwitch, 1));

%calculate average position and velocity for relative frame
qc = mean(qi0, 2); pc = mean(pi0, 2);

%sigma norm
sigma_norm = @(z,e) (1/e)*(sqrt(1 + e*norm(z)^2) - 1);
sigma_eps = @(z,e) z/sqrt(1+e*norm(z)^2);
d_alpha = sigma_norm(d, eps_parm);
r_alpha = sigma_norm(r, eps_parm);

%define phi_alpha(z);
h_phiAlpha = 0.2;
accumStep = 1e-3;
zVec = (0:0.1:25);
[phi_alpha, xi_alpha] = deal(zeros(length(zVec),1));
for ii = 1:length(zVec)
    z = zVec(ii);
    %integrate for xi
    if(z < d_alpha)
        ss = z:accumStep:d_alpha;
        sgn = -1.0;
    elseif z >= d_alpha
        ss = d_alpha:accumStep:z;
        sgn = 1.0;
    end
    tmp_xi_alpha = 0;
    for jj = 1:length(ss)
        s = ss(jj);
        phi_alpha(jj) = bump(s/r_alpha,h_phiAlpha) * action_base(a,b,c,s-d_alpha);
        tmp_xi_alpha = tmp_xi_alpha + phi_alpha(jj)*accumStep;
    end
    xi_alpha(ii) = sgn * tmp_xi_alpha;
end

if(0)
    figure('Name', 'Smooth Pairwise Potential');
    plot( zVec, xi_alpha, 'b' );
    xlabel("z"); ylabel("\psi_a");
end

%phi_beta(z)
h_phiBeta = 0.9;

%init member positions and velocities over time
[xi, vi] = deal(zeros(nDims,nNodes,nSamps));
xi(:,:,1) = qi0 - qc; vi(:,:,1) = pi0 - pc;

%translate gamma aganet into the relative frame
xd = qd - qc; vd = pd - pc;
xd_star = xd(:,1);

%propagate the dynamics
for tt = 1:nSamps
    
    %calculate the dynamics of the flocking alpha-agents
    for ii = 1:nNodes
        %check all neigbors for current node
        %Ni = {j el V : ||qj - qi|| < r}
        Ni = find(vecnorm(xi(:,:,tt) - xi(:, ii, tt)) < r);
        Ni = Ni(Ni~=ii); %remove self entry
        
        %if empty, don't control 
        if isempty(Ni)
            continue;
        end
        
        %if there are neighbors, form the control therms
        [fi_g, fi_d] = deal(zeros(nDims,1));
        for jj = Ni
            %gradient based term
            qjqi_signorm = sigma_norm(xi(:,jj,tt) - xi(:,ii,tt), eps_parm);
            phi_alpha = bump(qjqi_signorm/r_alpha,h_phiAlpha) * action_base(a,b,c,qjqi_signorm-d_alpha);
            fi_g = fi_g + phi_alpha * sigma_eps(xi(:,jj,tt) - xi(:,ii,tt), eps_parm);
            
            %velocity consensus term
            %calculate weight
            aij = bump(qjqi_signorm/r_alpha, h_phiAlpha);
            fi_d = fi_d + aij * (vi(:,jj,tt) - vi(:,ii,tt));
        end
        
        %check current switch time to adjust the navigation point
        if (tt*dT > tNavSwitch(switchIx))
            navIx = mod(switchIx,2)+1;
            xd_star = xd(:,navIx);
            if switchIx < length(tNavSwitch)
                switchIx = switchIx + 1;
            end
        end

        %create navigational feedback term
        fi_gamma = -c1 * (xi(:,ii,tt) - xd_star) - c2 * (vi(:,ii,tt) - vd);

        %calculate control signal
        uAlpha = fi_g + fi_d + fi_gamma;
        
        if(0) %rk4
            %integrate control into velocity
            vi(:,ii,tt+1) = rk4(vi(:,ii,tt), uAlpha, dT);

            %integrate velocity into position
            xi(:,ii,tt+1) = rk4(xi(:,ii,tt), vi(:,ii,tt+1), dT);

        elseif (1) %state space
            vi(:,ii,tt+1) = vi(:,ii,tt) + uAlpha * dT;
            xi(:,ii,tt+1) = xi(:,ii,tt) + vi(:,ii,tt+1)*dT;
        end

    end

    %plot agents
    if (tt == 1)
        figure('Name', 'Flocking Agents');
        xlabel("x [pos]"); ylabel("y [pos]");
    end
    if (mod(tt, 10) == 0)
        [deltaX_hat, deltaY_hat,deltaZ_hat,x0,y0,z0] = deal(zeros(nNodes,1));
        for ii = 1:nNodes
            % current position to current position + deltaT * v
            % x0(ii) = xi(1,ii,tt)+qc(1); x1 = x0(ii) + (vi(1,ii,tt)+pc(1)) * dT;
            % y0(ii) = xi(2,ii,tt)+qc(2); y1 = y0(ii) + (vi(2,ii,tt)+pc(2)) * dT;
            x0(ii) = xi(1,ii,tt); x1 = x0(ii) + vi(1,ii,tt) * dT;
            y0(ii) = xi(2,ii,tt); y1 = y0(ii) + vi(2,ii,tt) * dT;
            z0(ii) = xi(3,ii,tt); z1 = z0(ii) + vi(3,ii,tt) * dT;
            deltaX = x1-x0(ii); deltaY = y1-y0(ii); deltaZ = z1-z0(ii);
            deltaPos = [deltaX;deltaY;deltaZ]; normPos = norm(deltaPos);
            deltaX_hat(ii) = deltaX/normPos; deltaY_hat(ii) = deltaY/normPos;
            deltaZ_hat(ii) = deltaZ/normPos;
        end
        plot3(xd_star(1)+qc(1), xd_star(2)+qc(2), xd_star(3)+qc(3), 'xr', 'Linewidth', 3); hold on;
        quiver3(x0,y0,z0,deltaX_hat,deltaY_hat,deltaZ_hat,'color', 'b'); hold off;
        xlim([-50,250]); ylim([-10,60]); zlim([-10,110]);
        title(['t=', num2str(tt*dT), ' nav switch at t=', num2str(tNavSwitch(switchIx))]);
        drawnow; 
    end
    %calculate the dynamics of the gamma agents
    %for now I want the gamma agent to be static however in general the
    %agent can be a moving point

end

%% helper methods
function x = rk4(x0, f_x, h)
    
    %define integrator dynamics
    k1 = f_x;
    k2 = (f_x + h*k1/2);
    k3 = (f_x + h*k2/2);
    k4 = (f_x + h*k3);
    
    %approximate integral
    x = x0 + (h/6)*(k1 + 2*k2 + 2*k3 + k4);

end

function rho_hz = bump(z, h)
%implement the bump function from the flocking for multi agent paper

%check h
if (h < 0 || h > 1)
    error("Incorrect range for h entered, needs to be between 0 and 1");
end

%if h == 0, then bump is always zero
if (z < h && z >= 0)
    rho_hz = 1;
elseif (z >= h && z < 1)
    rho_hz = 0.5 * (1 + cos(pi*(z-h)/(1-h)));
else
    rho_hz = 0;
end

end

function phi_z = action_base(a,b,c,z)
    sig1 = @(z) z/sqrt(1 + z^2);
    phi_z = 0.5*((a+b)*sig1(z+c) + (a-b));
end