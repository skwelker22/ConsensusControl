%% consensus control flocking simulation
%startup & format
clear all; close all; home;

format long g;

%format
set(0, 'DefaultTextFontName', 'Arial', 'DefaultTextFontSize', 20, 'DefaultTextFontWEight', 'Bold');
set(0, 'DefaultAxesFontName', 'Arial', 'DefaultAxesFontSize', 20, 'DefaultAxesFontWeight', 'Bold', 'DefaultAxesLineWidth', 1.5);
set(0, 'DefaultLineLineWidth', 3, 'DefaultLineMarkerSize', 10);

%add libraries for controller design
addpath(genpath('/Users/swelker/Desktop/swelker/ConsensusControl'));

%% driver for controller
plotFlag=true;
dT = 1e-2;
truthDynamics=false;
nDims=2;

%create controller
[Tquad_ss2d, Kquad2d]=CreateQuadController(dT,nDims,truthDynamics,plotFlag);
