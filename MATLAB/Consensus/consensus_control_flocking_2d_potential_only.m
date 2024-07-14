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
eps_parm = 0.1; %fixed epsilon used for sigma-norm calculations
a = 5; %a parameter for potential function
b = a; %b parameter for potential function
c = abs(a-b)/sqrt(4*a*b); % c parameter for potential function
distVar = 7; %variance on the intial positions of the nodes
nNodes  = 10; %number of  agents in the networks
nDims   = 2; %number of dimensions of the problem
s = rng("default"); %fix random seed
qi0 = sqrt(distVar)*randn(nDims,nNodes); %initial agent positions
pRng = [-2,-1].^2; %range of initial agent velocities
pi0 = randi(fliplr(pRng), nDims, nNodes); %initial velocity
pi_angle = atan2d(pi0(2,:), pi0(1,:)); %relative heading angle

%define static gamma agent
qd = [[200; 30], [50;40]]; pd = [5;0];

%define switch times
nSwitch = 10; switchIx = 1;
tNavSwitch = unique(randi([10,tFinal], nSwitch, 1));

%calculate average position and velocity for relative frame
qc = mean(qi0, 2); pc = mean(pi0, 2);

%sigma norm
% sigma_norm = @(z,e) (1/e)*(sqrt(1 + e*norm(z)^2) - 1);
sigma_norm = @(z,e) (1/e)*(sqrt(1 + e.*vecnorm(z).^2) - 1);
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
    xlim([min(zVec), max(zVec)]);
end

%init member positions and velocities over time
[xi, vi] = deal(zeros(nDims,nNodes,nSamps));
xi(:,:,1) = qi0 - qc; vi(:,:,1) = pi0 - pc;

%define graph matrix to look at over time
qjqi_norm = zeros(nNodes, nNodes, nSamps);
EDGE_q = zeros(nNodes,1);
qiMinNorm = zeros(nNodes,nSamps);
uAlpha = zeros(nDims,nNodes,nSamps);
%propagate the dynamics
for tt = 1:nSamps
    
    EDGE_q_outer = 0; %outer sum for deviation energy
    
    %calculate the dynamics of the flocking alpha-agents
    for ii = 1:nNodes
        %check all neigbors for current node
        %Ni = {j el V : ||qj - qi|| < r}
        % Ni = find(vecnorm(xi(:,:,tt) - xi(:, ii, tt)) < r);
        Ni = find(sigma_norm(xi(:,:,tt) - xi(:, ii, tt), eps_parm) < r_alpha);
        Ni = Ni(Ni~=ii); %remove self entry
        EDGE_q(ii) = length(Ni);
        
        %if empty, don't control 
        if isempty(Ni)
            %if no control input, maintain velocity
            % xi(:,ii,tt+1) = xi(:,ii,tt) + vi(:,ii,tt).*dT;
            xi(:,ii,tt+1) = xi(:,ii,tt);
            vi(:,ii,tt+1) = vi(:,ii,tt);
            
            %calculate minimum norm with everyone
            minJNorm = 9999;
            for jj = 1:nNodes
                if (ii ~= jj)
                    qiqj_norm = norm(xi(:,jj,tt) - xi(:,ii,tt));
                    if(qiqj_norm < minJNorm)
                        minJNorm = qiqj_norm;
                    end
                end
            end
            qiMinNorm(ii,tt) = minJNorm;
            continue;
        end
        
        %if there are neighbors, form the control therms
        [fi_g, fi_d] = deal(zeros(nDims,1));
        EDGE_q_inner = 0; minJNorm = 9999;
        for jj = Ni
            %gradient based term
            qjqi_signorm = sigma_norm(xi(:,jj,tt) - xi(:,ii,tt), eps_parm);
            phi_alpha = bump(qjqi_signorm/r_alpha,h_phiAlpha) * action_base(a,b,c,qjqi_signorm-d_alpha);
            fi_g = fi_g + phi_alpha * sigma_eps(xi(:,jj,tt) - xi(:,ii,tt), eps_parm);

            %keep track of min norm
            qiqj_norm = norm(xi(:,jj,tt) - xi(:,ii,tt));
            if(qiqj_norm < minJNorm)
                minJNorm = qiqj_norm;
            end

            %velocity consensus term
            %calculate weight
            aij = bump(qjqi_signorm/r_alpha, h_phiAlpha);
            fi_d = fi_d + aij * (vi(:,jj,tt) - vi(:,ii,tt));

            %save off current graph matrix
            AGraph(ii,jj,tt) = aij;

            %inner sum term for deviation energy
            EDGE_q_inner = EDGE_q_inner + (norm(xi(:,jj,tt) - xi(:,ii,tt)) - d)^2;
        end % end jj
        
        %distace of the closest agent to evaluate the potential field
        %control
        qiMinNorm(ii,tt) = minJNorm;

        %calculate control signal
        disp(['Current u: ', num2str(fi_g'), ' at node: ', num2str(ii)]);
        uAlpha(:,ii,tt) = fi_g + fi_d;
        
        %state space breaks down when the dynamics are highly nonlinear,
        %but should work just fine for constant velocity assumption
        xi(:,ii,tt+1) = xi(:,ii,tt) + vi(:,ii,tt)*dT;
        vi(:,ii,tt+1) = vi(:,ii,tt) + uAlpha(:,ii,tt) * dT;

        %increment outer deviation energy sum
        EDGE_q_outer = EDGE_q_outer + EDGE_q_inner;

    end %end ii

    %calculate deviationi energyt to measure similarty to alpha-lattice
    EDGE_mag = sum(EDGE_q);
    E_q      = 1/(EDGE_mag + 1) * EDGE_q_outer;
    E_q_del  = EDGE_mag/(EDGE_mag + 1) * (eps_parm * d)^2;

    %plot agents
    if (tt == 1)
        f1 = figure('Name', 'Flocking Agents');
        xlabel("x [pos]"); ylabel("y [pos]");
        f3 = figure('Name', 'Min Distance from q_i to q_j');
        xlabel('i'); ylabel('||q_j-q_i||');
    end
    % if (mod(tt, 10) == 0)
        set(0, 'CurrentFigure', f1);
        [deltaX_hat, deltaY_hat,x0,y0] = deal(zeros(nNodes,1));
        for ii = 1:nNodes
            % current position to current position + deltaT * v
            x0(ii) = xi(1,ii,tt); x1 = x0(ii) + vi(1,ii,tt) * dT;
            y0(ii) = xi(2,ii,tt); y1 = y0(ii) + vi(2,ii,tt) * dT;
            deltaX = x1-x0(ii); deltaY = y1-y0(ii); deltaPos = [deltaX;deltaY];
            deltaX_hat(ii) = deltaX/norm(deltaPos); deltaY_hat(ii) = deltaY/norm(deltaPos);
            quiver(x0(ii),y0(ii),deltaX_hat(ii),deltaY_hat(ii),0,'color','b'); hold on;
            text(x0(ii), y0(ii), num2str(ii));
        end
        hold off;
        xlim([-50,50]); ylim([-30,30]);
        title(['t=', num2str(tt*dT), ...
               ' E(q): ', num2str(E_q), ...
               ' E(q)_\delta: ', num2str(E_q_del)]);
        drawnow; 
    % end
    if(mod(tt,.1/dT) == 0 && 1)
        %min distance
        set(0, 'CurrentFigure', f3);
        plot(qiMinNorm(:,tt)); yline(r, 'LineWidth',4);
        drawnow;
    end

    %plot controls and states
    if(0)
        ttPlot = 1:tt; uCCC = winter(nNodes);
        %controls
        ax = figure('Name', 'Controls vs Time'); legCell = cell(nNodes,1);
        for ii = 1:nNodes
            subplot(211);
            plot( ttPlot, squeeze(uAlpha(1,ii,ttPlot)), 'Color', uCCC(ii,:)); hold on;
            subplot(212);
            plot( ttPlot, squeeze(uAlpha(2,ii,ttPlot)), 'Color', uCCC(ii,:)); hold on;
            legCell{ii} = num2str(ii);
        end
        hold off; xlabel('Time [sec]');
        subplot(211); ylabel("u_x"); legend(legCell);
        subplot(212); ylabel("u_y");
        %states
        figure('Name', 'States vs. Time');
        for ii = 1:nNodes
            subplot(221);
            plot( ttPlot, squeeze(xi(1,ii,ttPlot)), 'Color', uCCC(ii,:)); hold on;
            subplot(222);
            plot( ttPlot, squeeze(xi(2,ii,ttPlot)), 'Color', uCCC(ii,:)); hold on;
            subplot(223);
            plot( ttPlot, squeeze(vi(1,ii,ttPlot)), 'Color', uCCC(ii,:)); hold on;
            subplot(224);
            plot( ttPlot, squeeze(vi(2,ii,ttPlot)), 'Color', uCCC(ii,:)); hold on;
        end
        hold off; xlabel('Time [sec]');
        subplot(221); ylabel("p_x"); subplot(222); ylabel("p_y");
        subplot(223); ylabel("v_x"); subplot(224); ylabel("v_y");
    end

end %end tt

%% position/velocity and control plots
figure('Name','Control');
ttplt = (1:tt-1);
ccc = turbo(nNodes);
for nn = 1:nNodes
    subplot(211);
    plot(dT.*ttplt, squeeze(uAlpha(1,nn,1:tt-1)), 'Color', ccc(nn,:)); hold on;
    subplot(212);
    plot(dT.*ttplt, squeeze(uAlpha(2,nn,1:tt-1)), 'Color', ccc(nn,:)); hold on;
end
subplot(211); ylabel("u_{\alpha,x}");
subplot(212); xlabel("Time [sec]"); ylabel("u_{\alpha,y}");

% figure('Name','Positions');
% for nn = 1:nNodes
%     subplot(211);
%     plot(dT.*ttplt, squeeze(xi(1,nn,1:tt-1)), 'Color', ccc(nn,:)); hold on;
%     subplot(212);
%     plot(dT.*ttplt, squeeze(xi(2,nn,1:tt-1)), 'Color', ccc(nn,:)); hold on;
% end
% subplot(211); ylabel("pos_x");
% subplot(212); xlabel("Time [sec]"); ylabel("pos_y");
% 
% figure('Name','Velocities');
% for nn = 1:nNodes
%     subplot(211);
%     plot(dT.*ttplt, squeeze(vi(1,nn,1:tt-1)), 'Color', ccc(nn,:)); hold on;
%     subplot(212);
%     plot(dT.*ttplt, squeeze(vi(2,nn,1:tt-1)), 'Color', ccc(nn,:)); hold on;
% end
% subplot(211); ylabel("vel_x");
% subplot(212); xlabel("Time [sec]"); ylabel("vel_y");

%% helper methods
function rho_hz = bump(z, h)
%implement the bump function from the flocking for multi agent paper

%check h
if (h < 0 || h > 1)
    error("Incorrect range for h entered, needs to be between 0 and 1");
end

%if h == 0, then bump is always zero
if (z < h && z >= 0)
    rho_hz = 1;
elseif (z >= h && z <= 1)
    rho_hz = 0.5 * (1 + cos(pi*(z-h)/(1-h)));
else
    rho_hz = 0;
end

end

function phi_z = action_base(a,b,c,z)
    sig1 = @(z) z/sqrt(1 + z^2);
    phi_z = 0.5*((a+b)*sig1(z+c) + (a-b));
end


