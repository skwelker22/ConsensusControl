function [sysMag, sysPhase, sysOmega] = getBodeData(sys)
%gets the bode magnitude, phase and angle data

%and then squeezes out the 3rd dimension for use in computations
[tmpMag, tmpPhase, tmpOmega] = bode(sys);

%set tmp to output
sysMag   = squeeze(tmpMag(1,:,:)); 
sysPhase = squeeze(tmpPhase(1,:,:));
sysOmega = tmpOmega;

end
