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
% nNodes  = 150; %number of  agents in the networks
distVar = 150;
nNodes  = 20;
nDims   = 2; %number of dimensions of the problem
% s = rng("default"); %fix random seed
s=rng(3);
qi0 = sqrt(distVar)*randn(nDims,nNodes); %initial agent positions
pRng = [-2,-1].^2; %range of initial agent velocities
pi0 = randi(fliplr(pRng), nDims, nNodes); %initial velocity
pi_angle = atan2d(pi0(2,:), pi0(1,:)); %relative heading angle

%calculate average position and velocity for relative frame
qc0 = mean(qi0, 2); pc0 = mean(pi0, 2);

%define static gamma agent
% qd = [[200; 30], [50;40]]; pd = [5;0];
qd = [[0;0], [200;30]]; pd = [5;0];

%define switch times
nSwitch = 10; switchIx = 1;
% tNavSwitch = unique(randi([10,tFinal], nSwitch, 1));
tNavSwitch = [50, 55];

%sigma norm
% sigma_norm = @(z,e) (1/e)*(sqrt(1 + e*norm(z)^2) - 1);
sigma_norm = @(z,e) (1/e)*(sqrt(1 + e.*vecnorm(z).^2) - 1);
sigma_eps = @(z,e) z/sqrt(1+e*norm(z)^2);
d_alpha = sigma_norm(d, eps_parm);
r_alpha = sigma_norm(r, eps_parm);

%define phi_alpha(z);
h_phiAlpha = 0.2;
accumStep = 1e-3;
zStep = 0.1; zVec = (0:zStep:25);
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
    plot( zVec, xi_alpha, 'b' ); hold on;
    plot( zVec, -[diff(xi_alpha);inf]./zStep, 'r');
    xlabel("z"); ylabel("\psi_a"); legend('V(z)', '-\nablaV(z)');
    xlim([min(zVec), max(zVec)]);
end

%phi_beta(z)
h_phiBeta = 0.9;

%init member positions and velocities over time
[xi, vi] = deal(zeros(nDims,nNodes,nSamps));
[qc, pc] = deal(zeros(nDims, nSamps));
xi(:,:,1) = qi0 - qc0; vi(:,:,1) = pi0 - pc0;
qc(:,1) = qc0; pc(:,1) = pc0;

%translate gamma agent into the relative frame
xd = qd - qc0; vd = pd - pc0;
xd_star = xd(:,1);

%define graph matrix to look at over time
[AGraph,qjqi_norm] = deal(zeros(nNodes, nNodes, nSamps));
EDGE_q = zeros(nNodes,1);
qiMinNorm = zeros(nNodes,nSamps);
uAlpha = zeros(nDims,nNodes,nSamps);
[fi_g_plt, fi_d_plt, fi_gamma_plt] = deal(zeros(nDims,nNodes,nSamps));
xd_plot = zeros(nDims, nSamps); E_q_bar = zeros(nSamps,1);
%propagate the dynamics
for tt = 1:nSamps
    %check current switch time to adjust the navigation point
    if (tt*dT > tNavSwitch(switchIx))
        navIx = mod(switchIx,2)+1;
        xd_star = xd(:,navIx);
        if switchIx < length(tNavSwitch)
            switchIx = switchIx + 1;
        end
    end
    xd_plot(:,tt) = xd_star;

    %initialize outer sum for deviation energy
    EDGE_q_outer = 0; 
    %calculate the dynamics of the flocking alpha-agents
    for ii = 1:nNodes
        %init min norm
        minJNorm = 9999;

        %check all neigbors for current node
        %Ni = {j el V : ||qj - qi||_sigma < r_sigma}
        Ni = find(sigma_norm(xi(:,:,tt) - xi(:, ii, tt), eps_parm) < r_alpha);
        Ni = Ni(Ni~=ii); %remove self entry
        EDGE_q(ii) = length(Ni);
        
        %if empty, only apply the navigational term
        if isempty(Ni)
            %create navigational feedback term
            %this should apply even when no connections are available as
            %this would programmed into the mission
            fi_gamma = -c1 * (xi(:,ii,tt) - xd_star) - c2 * (vi(:,ii,tt) - vd);
            % fi_gamma = zeros(2,1);
            uAlpha(:,ii,tt) = fi_gamma;

            %if no control input, maintain velocity
            xi(:,ii,tt+1) = xi(:,ii,tt) + vi(:,ii,tt) * dT;
            vi(:,ii,tt+1) = vi(:,ii,tt) + uAlpha(:,ii,tt) * dT;

            %calculate minimum norm with everyone
            for jj = 1:nNodes
                if (ii ~= jj)
                    qiqj_norm = sigma_norm(xi(:,jj,tt) - xi(:,ii,tt), eps_parm);
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
        EDGE_q_inner = 0; 
        for jj = Ni
            %gradient based term
            qjqi_signorm = sigma_norm(xi(:,jj,tt) - xi(:,ii,tt), eps_parm);
            phi_alpha = bump(qjqi_signorm/r_alpha,h_phiAlpha) * action_base(a,b,c,qjqi_signorm-d_alpha);
            fi_g = fi_g + phi_alpha * sigma_eps(xi(:,jj,tt) - xi(:,ii,tt), eps_parm);

            %keep track of min norm
            qiqj_norm = sigma_norm(xi(:,jj,tt) - xi(:,ii,tt), eps_parm);
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
            EDGE_q_inner = EDGE_q_inner + ...
                (sigma_norm(xi(:,jj,tt) - xi(:,ii,tt),eps_parm) - d_alpha)^2;
        end % end jj
        
        %keep track of the min norm for the closest agent to track the
        %structure of the alpha lattice
        qiMinNorm(ii,tt) = minJNorm;
        
        %create navigational feedback term
        fi_gamma = -c1 * (xi(:,ii,tt) - xd_star) - c2 * (vi(:,ii,tt) - vd);
        % fi_gamma = zeros(2,1);

        %calculate control signal
        disp(['Current u: ', num2str(fi_g'), ' at node: ', num2str(ii)]);
        uAlpha(:,ii,tt) = fi_g + fi_d + fi_gamma;
        fi_g_plt(:,ii,tt) = fi_g;
        fi_d_plt(:,ii,tt) = fi_d;
        fi_gamma_plt(:,ii,tt) = fi_gamma;
        
        %two different implementations of the dynamics
        %rk4 integration and state space
        %state space breaks down when the dynamics are highly nonlinear,
        %but should work just fine for constant velocity assumption
        if(0) %rk4
            %integrate control into velocity
            vi(:,ii,tt+1) = rk4(vi(:,ii,tt), uAlpha, dT);

            %integrate velocity into position
            xi(:,ii,tt+1) = rk4(xi(:,ii,tt), vi(:,ii,tt+1), dT);

        elseif (1) %state space
            xi(:,ii,tt+1) = xi(:,ii,tt) + vi(:,ii,tt) * dT;
            vi(:,ii,tt+1) = vi(:,ii,tt) + uAlpha(:,ii,tt) * dT;
        end

        %increment outer deviation energy sum
        EDGE_q_outer = EDGE_q_outer + EDGE_q_inner;

    end %end ii

    %update center of mass
    qc(:,tt+1) = qc(:,tt) + pc(:,tt) * dT;
    pc(:,tt+1) = pc(:,tt) + mean(squeeze(uAlpha(:,:,tt)),2) * dT;

    %calculate deviation energy to measure similarty to alpha-lattice
    EDGE_mag = sum(EDGE_q);
    E_q      = 1/(EDGE_mag + 1) * EDGE_q_outer;
    E_q_bar(tt) = E_q/d_alpha^2;
    % E_q_del  = EDGE_mag/(EDGE_mag + 1) * (1e-1 * d_alpha)^2;
    E_q_del = (1e-1 * d_alpha)^2;

    %plot agents
    if (tt == 1)
        f1 = figure('Name', 'Flocking Agents');
        xlabel("x [pos]"); ylabel("y [pos]");
        % f2 = figure('Name', 'Graph Adjacency Matrix');
        % xlabel("j"); ylabel("i");
        f3 = figure('Name', 'Min Distance from q_i to q_j');
        xlabel('i'); ylabel('||q_j-q_i||');
    end
    if (mod(tt, 10) == 0)
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
        scatter(xd_star(1), xd_star(2), 'xr', 'Linewidth', 3); hold off
        % xlim([-50,250]); ylim([-10,60]);
        xlim([-50,50]); ylim([-50,50]);
        title(['t=', num2str(tt*dT), ...
               ' nav switch at t=', num2str(tNavSwitch(switchIx)), ...
               ' E(q): ', num2str(E_q), ...
               ' E(q)_\delta: ', num2str(E_q_del)]);
        drawnow; 
    end
    if(mod(tt,.1/dT) == 0 && 1)
        %graph adjaceny matrix
        % set(0, 'CurrentFigure', f2); 
        % allNodes = 1:nNodes; [XX,YY] = meshgrid(allNodes, allNodes);
        % h = surf(XX, YY, AGraph(:,:,tt)); 
        % colorbar; view(0,90);
        % drawnow;
        %min distance
        set(0, 'CurrentFigure', f3);
        plot(qiMinNorm(:,tt), 'b'); 
        yline(r, 'LineWidth',4);
        yline(d_alpha, 'LineWidth', 4, 'Color', [0.49,0.18,0.56]);
        ylim([0, 2*r]); legend('q_{ij, min}', 'r', 'd_\alpha');
        drawnow;
    end

    %plot controls and states
    if(0)
        ttIx = 1:tt; ttPlot = dT.*(ttIx-1);
        switchPlot = tNavSwitch(tNavSwitch<=ttPlot(end));
        if isempty(switchPlot), switchPlot = 0; end
        uCCC = winter(nNodes);
        nodeSubset = [5,15];
        %controls
        ax = figure('Name', 'Controls vs Time'); legCell = cell(nNodes,1);
        hh = zeros(nNodes,1);
        for ii = 1:nNodes
            subplot(211);
            hh(ii) = plot( ttPlot, squeeze(uAlpha(1,ii,ttIx)), 'Color', uCCC(ii,:)); hold on;
            subplot(212);
            plot( ttPlot, squeeze(uAlpha(2,ii,ttIx)), 'Color', uCCC(ii,:)); hold on;
            legCell{ii} = num2str(ii);
        end
        subplot(211); xline(switchPlot, 'k');
        hold off; xlabel('Time [sec]'); ylabel("u_x"); 
        lgd = legend(hh,legCell); fontsize(lgd, 9, 'points');
        subplot(212); ylabel("u_y"); xline(switchPlot, 'k');
        %states
        figure('Name', 'States vs. Time');
        for ii = 1:nNodes
            subplot(221);
            plot( ttPlot, squeeze(xi(1,ii,ttIx)), 'Color', uCCC(ii,:)); hold on;
            subplot(222);
            plot( ttPlot, squeeze(xi(2,ii,ttIx)), 'Color', uCCC(ii,:)); hold on;
            subplot(223);
            plot( ttPlot, squeeze(vi(1,ii,ttIx)), 'Color', uCCC(ii,:)); hold on;
            subplot(224);
            plot( ttPlot, squeeze(vi(2,ii,ttIx)), 'Color', uCCC(ii,:)); hold on;
        end
        %plot navigational point
        subplot(221); plot( ttPlot, xd_plot(1,ttIx), 'k' );
        subplot(222); plot( ttPlot, xd_plot(2,ttIx), 'k' );
        hold off; xlabel('Time [sec]');
        subplot(221); ylabel("p_x"); subplot(222); ylabel("p_y");
        subplot(223); ylabel("v_x"); subplot(224); ylabel("v_y");
        %plots for node subset only for debug
        %subset controls
        figure('Name', 'Subset Controls vs Time'); legCell2 = cell(length(nodeSubset),1);
        hhs = zeros(length(nodeSubset),1); nsIx = 1;
        for ss = nodeSubset
            subplot(211);
            hhs(nsIx) = plot( ttPlot, squeeze(uAlpha(1,ss,ttIx)), 'Color', uCCC(ss,:)); hold on;
            subplot(212);
            plot( ttPlot, squeeze(uAlpha(2,ss,ttIx)), 'Color', uCCC(ss,:)); hold on;
            legCell2{nsIx} = num2str(ss);
            nsIx = nsIx + 1;
        end
        subplot(211); xline(switchPlot, 'k');
        hold off; xlabel('Time [sec]'); ylabel("u_x"); title('Node Subset u');
        lgd = legend(hhs,legCell2); fontsize(lgd, 9, 'points');
        subplot(212); ylabel("u_y"); xline(switchPlot, 'k');
        %subset with control components
        figure('Name', ['Subset Control Components vs Time: Nodes: ',num2str(nodeSubset)]);
        legCell2 = cell(length(nodeSubset),1);
        hhs = zeros(length(nodeSubset),1); nsIx = 1;
        for ss = nodeSubset
        % for ss = 1:nNodes
        % for ss = 1:2
            subplot(321);
            hhs(nsIx) = plot( ttPlot, squeeze(fi_g_plt(1,ss,ttIx)), 'Color', uCCC(ss,:)); hold on;
            subplot(323);
            plot( ttPlot, squeeze(fi_d_plt(1,ss,ttIx)), 'Color', uCCC(ss,:)); hold on;
            subplot(325);
            plot( ttPlot, squeeze(fi_gamma_plt(1,ss,ttIx)), 'Color', uCCC(ss,:)); hold on;
            subplot(322);
            plot( ttPlot, squeeze(fi_g_plt(2,ss,ttIx)), 'Color', uCCC(ss,:)); hold on;
            subplot(324);
            plot( ttPlot, squeeze(fi_d_plt(2,ss,ttIx)), 'Color', uCCC(ss,:)); hold on;
            subplot(326);
            plot( ttPlot, squeeze(fi_gamma_plt(2,ss,ttIx)), 'Color', uCCC(ss,:)); hold on;
            legCell2{nsIx} = num2str(ss);
            nsIx = nsIx + 1;
        end
        subplot(321); xline(switchPlot, 'k');
        hold off; xlabel('Time [sec]'); ylabel("f_{g,x}");
        lgd = legend(hhs,legCell2); fontsize(lgd, 9, 'points');
        subplot(323); ylabel("f_{d,x}"); xline(switchPlot, 'k');
        subplot(325); ylabel("f_{\gamma,x}"); xline(switchPlot, 'k');
        subplot(322); ylabel("f_{g,y}"); xline(switchPlot, 'k');
        subplot(324); ylabel("f_{d,y}"); xline(switchPlot, 'k');
        subplot(326); ylabel("f_{\gamma,y}"); xline(switchPlot, 'k');
        %subset states
        figure('Name', 'Subset States vs. Time');
        for ss = nodeSubset
            subplot(221);
            plot( ttPlot, squeeze(xi(1,ss,ttIx)), 'Color', uCCC(ss,:)); hold on;
            subplot(222);
            plot( ttPlot, squeeze(xi(2,ss,ttIx)), 'Color', uCCC(ss,:)); hold on;
            subplot(223);
            plot( ttPlot, squeeze(vi(1,ss,ttIx)), 'Color', uCCC(ss,:)); hold on;
            subplot(224);
            plot( ttPlot, squeeze(vi(2,ss,ttIx)), 'Color', uCCC(ss,:)); hold on;
        end
        %plot navigational point
        subplot(221); plot( ttPlot, xd_plot(1,ttIx), 'k' );
        subplot(222); plot( ttPlot, xd_plot(2,ttIx), 'k' );
        hold off; xlabel('Time [sec]');
        subplot(221); ylabel("p_x"); subplot(222); ylabel("p_y");
        subplot(223); ylabel("v_x"); subplot(224); ylabel("v_y");

        %sig norm of the subset
        zzPlt = squeeze(xi(:,nodeSubset(1),ttIx))-squeeze(xi(:,nodeSubset(2),ttIx));
        zzSigPlt = sigma_norm(zzPlt, eps_parm);
        figure('Name', 'Difference vs. Sigma');
        plot( zzPlt(1,:), zzSigPlt, 'Color', uCCC(1,:) ); hold on;
        plot( zzPlt(2,:), zzSigPlt, 'Color', uCCC(end,:) );
        xlabel("x_i - x_j"); ylabel('||x_i - x_j||_\sigma');

        figure('Name', 'Sigma vs. Time & Diff vs. Time');
        plot( ttPlot, zzSigPlt, 'k'); hold on;
        plot( ttPlot, zzPlt(1,:), 'Color', uCCC(1,:));
        plot( ttPlot, zzPlt(2,:), 'Color', uCCC(end,:));
        xlabel("Time [sec]"); ylabel("||x_i - x_j||_\sigma");
        legend('||x_5 - x_{15}||_\sigma', 'x(1)_5 - x(1)_{15}', 'x(2)_5 - x(2)_{15}');
        
        %alpha-flock verification plots
        connectivity = zeros(tt,1); cohesion = zeros(tt,1);
        K_v_bar = zeros(tt,1);
        for nn = 1:tt
            %connectivity
            L_graph = diag(sum(AGraph(:,:,nn),2)) - AGraph(:,:,nn);
            connectivity(nn) = 1/(nNodes-1) * rank(L_graph);
            %cohesion
            cohesion(nn) = max(vecnorm(xi(:,:,nn)));
            %normalized velocity mismatch
            K_v_bar(nn) = (0.5/nNodes) * sum(vecnorm(vi(:,:,nn)).^2);
        end
        endLim = 1000;
        figure('Name', '\alpha-flock verification');
        subplot(221); plot( ttIx, connectivity, 'b'); ylabel("connectivity");
        xlim([0,endLim]); ylim([0,1]);
        subplot(222); plot( ttIx, cohesion, 'b' ); ylabel("cohesion radius");
        xlim([0,endLim]);
        subplot(223); plot( ttIx, E_q_bar(ttIx), 'b'); ylabel("deviation energy");
        xlim([0,endLim]);
        subplot(224); plot( ttIx, K_v_bar, 'b'); ylabel("velocity mismatch");
        xlim([0,endLim]);
    end

end %end tt

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

