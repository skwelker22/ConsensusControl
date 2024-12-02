function [ccfg] = GetConsensusConfig(rndSeed)

%set random seed
s=rng(rndSeed);

%set configurations
ccfg.dT      = 0.001;
ccfg.tFinal  = 15; %seconds
ccfg.nSamps  = ccfg.tFinal/ccfg.dT;
ccfg.distVar = 1; %this works
ccfg.nNodes  = 1;
ccfg.nDims   = 2; %number of dimensions of the problem

end