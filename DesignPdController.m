%% design PD controller
function C_PD = DesignPdController(Ps, omegaC, PM)
%function to design PD controller
% Ps : plant dynamics as a transfer function
% omegaC : crossover frequency ~ BW/1.5 - rad/sec
% PM : targeted phase margin for the design - deg

%form of the controller
%C = k (s + z)/(tau*s+1)
tau = 1/(10*omegaC);

%solve for derivative zero using the phase margin
[magP, phaseP, omegaP] = getBodeData(Ps);

%find omega index closest to omegaC
omegaCIx = find(omegaC>=omegaP, 1, 'last');

%phase of the plant at this index
angleP_deg = phaseP(omegaCIx);

%phase of the derivative filter
angleTau_deg = atan2d(tau*omegaC,1);

%calculate z
z = tand(-180 + PM - angleP_deg + angleTau_deg)/omegaC;

%calculate k
%find mag p closest to phase margin
magP_omegaC = magP(omegaCIx);

k = sqrt((tau*omegaC)^2+1)/sqrt(omegaC^2 + z) * (1/magP_omegaC);

%PD controller
C_PD = k*tf([1, z], [tau, 1]);

end