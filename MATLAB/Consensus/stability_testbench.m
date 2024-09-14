%%
eps_parm = 0.1; %fixed epsilon used for sigma-norm calculations
c1=1; c2=1;
distVar = 2500; %variance on the intial positions of the nodes
nNodes  = 150; %number of  agents in the networks
% distVar = 150;
% nNodes = 20;
% nNodes  = 30;
nDims   = 2; %number of dimensions of the problem
% s=rng(3);
qi0 = sqrt(distVar)*randn(nDims,nNodes); %initial agent positions
pRngX = [-2,2]; pRngY = [-1,1];
pi0X = randi(pRngX, 1, nNodes);
pi0Y = randi(pRngY, 1, nNodes);
pi0  = [pi0X; pi0Y];

%calculate average position and velocity for relative frame
qc0 = mean(qi0, 2); pc0 = mean(pi0, 2);

%display c* which is the initial energy in the system
cstar=CalculatePotential(0,eps_parm);
disp(["c*="+num2str(cstar)]);

%calculate hamiltonian
%init member positions and velocities over time
[xi, vi] = deal(zeros(nDims,nNodes));
xi(:,:) = qi0 - qc0; vi(:,:) = pi0 - pc0;

%calculate the dynamics of the flocking alpha-agents
allNodes=(1:nNodes)'; Vx=0;
for ii = 1:nNodes
    qjVec=allNodes(allNodes~=ii)';
    Vx_inner=0;
    for jj = qjVec
        qjqi=sigma_norm(xi(:,jj)-xi(:,ii),eps_parm);
        Vx_inner=Vx_inner+CalculatePotential(qjqi,eps_parm);
    end
    Vx=Vx+Vx_inner;
end
Vx=0.5*Vx;

%calculate J(x),K(v)
Jx=c1/2*sum(vecnorm(xi).^2);
Kv=0.5*sum(vecnorm(vi).^2);
Hxv0=Vx+Jx+Kv;
disp("H(x(0),v(0))="+num2str(Hxv0));

%calculate k
khat=ceil(Hxv0/cstar-1);
disp("k="+num2str(khat));

function xi_alpha=CalculatePotential(z, eps_parm)
    h_phiAlpha = 0.2;
    accumStep = 1e-2;
    d = 7; %constraint distance agent to agent
    r = 1.2 * d; %ball radius which defines the max distance for the proximity net
    d_alpha = sigma_norm(d, eps_parm);
    r_alpha = sigma_norm(r, eps_parm);
    % a = 10; b = 20;
    a = 5; b = 5;
    % a = 1; b = 1;
    c = abs(a-b)/sqrt(4*a*b); % c parameter for potential function
    
    %integrate for xi
    if(z < d_alpha)
        ss = z:accumStep:d_alpha;
        sgn = -1.0;
    elseif z >= d_alpha
        ss = d_alpha:accumStep:z;
        sgn = 1.0;
    end
    % tmp_xi_alpha = 0;
    % for jj = 1:length(ss)
    %     s = ss(jj);
    %     phi_alpha = bump(s/r_alpha,h_phiAlpha) * action_base(a,b,c,s-d_alpha);
    %     tmp_xi_alpha = tmp_xi_alpha + phi_alpha*accumStep;
    % end
    % xi_alpha = sgn * tmp_xi_alpha;
    phi_alpha=@(x)(sgn*bump(x/r_alpha,h_phiAlpha)*action_base(a,b,c,x-d_alpha));
    xi_alpha=integral(phi_alpha,ss(1),ss(end),'ArrayValued',true);
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

function sig_norm=sigma_norm(z,e)
    % sigma_norm = @(z,e) (1/e)*(sqrt(1 + e.*vecnorm(z).^2) - 1);
    sig_norm=(1/e)*(sqrt(1 + e.*vecnorm(z).^2) - 1);
end